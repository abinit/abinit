!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pawdij
!! NAME
!!  m_pawdij
!!
!! FUNCTION
!!  This module contains several routines used to compute the PAW pseudopotential
!!  strengths Dij. The Dijs define the non-local PAW operator:
!!         VNL = Sum_ij [ Dij |pi><pj| ],  with pi, pj= projectors
!!
!! COPYRIGHT
!! Copyright (C) 2013-2019 ABINIT group (MT, FJ, BA, JWZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

MODULE m_pawdij

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING

 use m_paral_atom,   only : get_my_atmtab, free_my_atmtab
 use m_paw_io,       only : pawio_print_ij
 use m_pawang,       only : pawang_type
 use m_pawrad,       only : pawrad_type, pawrad_deducer0, simp_gen, nderiv_gen
 use m_pawtab,       only : pawtab_type
 use m_paw_an,       only : paw_an_type
 use m_paw_ij,       only : paw_ij_type, paw_ij_print
 use m_pawfgrtab,    only : pawfgrtab_type
 use m_pawrhoij,     only : pawrhoij_type
 use m_paw_finegrid, only : pawgylm, pawexpiqr
 use m_paw_sphharm,  only : slxyzs

 implicit none

 private

!public procedures.
 public :: pawdij           ! Dij total
 public :: pawdijhartree    ! Dij Hartree
 public :: pawdijfock       ! Dij Fock exact-exchange
 public :: pawdijxc         ! Dij eXchange-Correlation (using (r,theta,phi) grid)
 public :: pawdijxcm        ! Dij eXchange-Correlation (using (l,m) moments)
 public :: pawdijhat        ! Dij^hat (compensation charge contribution)
 public :: pawdijnd         ! Dij nuclear dipole
 public :: pawdijso         ! Dij spin-orbit
 public :: pawdiju          ! Dij LDA+U
 public :: pawdiju_euijkl   ! Dij LDA+U, using pawrhoij instead of occupancies
 public :: pawdijexxc       ! Dij local exact-exchange
 public :: pawdijfr         ! 1st-order frozen Dij
 public :: pawpupot         ! On-site LDA+U potential
 public :: pawxpot          ! On-site local exact-exchange potential
 public :: symdij           ! Symmetrize total Dij or one part of it
 public :: symdij_all       ! Symmetrize all contributions to Dij
 public :: pawdij_gather    ! Perform a allgather operation on Dij
 public :: pawdij_print_dij ! Print out a Dij matrix
!!***

CONTAINS

!===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawdij
!! NAME
!! pawdij
!!
!! FUNCTION
!! Compute the pseudopotential strengths Dij of the PAW non local operator as sum of
!! several contributions. Can compute first-order strenghts Dij for RF calculations.
!! This routine is a driver calling, for each contribution to Dij, a specific
!! routines.
!! Within standard PAW formalism, Dij can be decomposd as follows:
!!      Dij = Dij_atomic + Dij_Hartree + Dij_XC + Dij^hat
!! In case of additional approximations, several other terms can appear:
!!      Dij_LDA+U, Dij_spin-orbit, Dij_local-exact-exchange, Dij_Fock...
!!
!! INPUTS
!!  cplex=1 if no phase is applied (GS), 2 if a exp(-iqr) phase is applied (Response Function at q<>0)
!!  enunit=choice for units of output Dij
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!  [hyb_mixing, hyb_mixing_sr]= -- optional-- mixing factors for the global (resp. screened) XC hybrid functional
!!  ipert=index of perturbation (used only for RF calculation ; set ipert<=0 for GS calculations.
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  nfft=number of real space grid points (for current proc)
!!  nfftot=total number of real space grid points
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms in unit cell.
!!  paw_an(my_natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  paw_ij(my_natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(my_natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawspnorb=flag: 1 if spin-orbit coupling is activated
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  qphon(3)=wavevector of the phonon
!!  spnorbscl=scaling factor for spin-orbit coupling
!!  ucvol=unit cell volume
!!  vtrial(cplex*nfft,nspden)=GS potential on real space grid
!!  vxc(cplex*nfft,nspden)=XC potential (Hartree) on real space grid
!!  xred(3,my_natom)= reduced atomic coordinates
!!  ======== Optional arguments ==============
!!  Parallelism over atomic sites:
!!    mpi_atmtab(:)=indexes of the atoms treated by current proc
!!    comm_atom=MPI communicator over atoms
!!    mpi_comm_grid=MPI communicator over real space grid points
!!  Application of a potential energy shift on atomic sites:
!!    natvshift=number of atomic potential energy shifts (per atom) ; default=0
!!    atvshift(natvshift,nsppol,natom)=potential energy shift for lm channel & spin & atom
!!    fatvshift=factor that multiplies atvshift
!!  Electrons-positron 2-component DFT:
!!    electronpositron_calctype=type of calculation for electron-positron 2component-DFT:
!!       0: standard DFT (no positron) ; default value
!!       1: positron  in the constant electrons potential
!!       2: electrons in the constant positron potential
!!    electronpositron_pawrhoij(my_natom) <type(pawrhoij_type)>=
!!       PAW occupation matrix of the "constant" particle(s)
!!       (electrons if calctype=1, positron if calctype=2)
!!    electronpositron_lmselect(lmmax,my_natom)=
!!       Flags selecting the non-zero LM-moments of on-site densities
!!       for the "constant" particle(s)
!!       (electrons if calctype=1, positron if calctype=2)
!!
!! OUTPUT
!!  paw_ij(iatom)%dij(cplex_dij*qphase*lmn2_size,ndij)= total Dij terms (GS calculation, ipert=0)
!!                                                   total 1st-order Dij terms (RF ccalc., ipert>0)
!!  May be complex if cplex_dij=2
!!        dij(:,1) contains Dij^up-up
!!        dij(:,2) contains Dij^dn-dn
!!        dij(:,3) contains Dij^up-dn (only if nspinor=2)
!!        dij(:,4) contains Dij^dn-up (only if nspinor=2)
!!  May also compute paw_ij(iatom)%dij0,paw_ij(iatom)%dijhartree,paw_ij(iatom)%dijxc,
!!                   paw_ij(iatom)%dijxc_hat,paw_ij(iatom)%dijxc_val,
!!                   paw_ij(iatom)%dijhat,paw_ij(iatom)dijso,
!!                   paw_ij(iatom)%dijU,paw_ij(iatom)%dijexxc,paw_ij(iatom)%dijfock
!!
!! NOTES
!!  Response function calculations:
!!    In order to compute first-order Dij, paw_an (resp. paw_ij) datastructures
!!    must contain first-order quantities, namely paw_an1 (resp. paw_ij1).
!!
!! PARENTS
!!      bethe_salpeter,dfpt_scfcv,respfn,scfcv,screening,sigma
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdij(cplex,enunit,gprimd,ipert,my_natom,natom,nfft,nfftot,nspden,ntypat,&
&          paw_an,paw_ij,pawang,pawfgrtab,pawprtvol,pawrad,pawrhoij,pawspnorb,pawtab,&
&          pawxcdev,qphon,spnorbscl,ucvol,charge,vtrial,vxc,xred,&
&          electronpositron_calctype,electronpositron_pawrhoij,electronpositron_lmselect,&
&          atvshift,fatvshift,natvshift,nucdipmom,&
&          mpi_atmtab,comm_atom,mpi_comm_grid,hyb_mixing,hyb_mixing_sr,mgga)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,enunit,ipert,my_natom,natom,nfft,nfftot
 integer,intent(in) :: nspden,ntypat,pawprtvol,pawspnorb,pawxcdev
 integer,optional,intent(in) :: electronpositron_calctype
 integer,optional,intent(in) :: comm_atom,mpi_comm_grid,natvshift
 integer,optional,intent(in) :: mgga
 real(dp),intent(in) :: spnorbscl,ucvol,charge
 real(dp),intent(in),optional ::fatvshift,hyb_mixing,hyb_mixing_sr
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 logical,optional,intent(in) :: electronpositron_lmselect(:,:)
 real(dp),intent(in) :: gprimd(3,3),qphon(3)
 real(dp),intent(in) ::  vxc(:,:),xred(3,natom)
 real(dp),intent(in),target :: vtrial(cplex*nfft,nspden)
 real(dp),intent(in),optional :: atvshift(:,:,:)
 real(dp),intent(in),optional :: nucdipmom(3,my_natom)
 type(paw_an_type),intent(in) :: paw_an(my_natom)
 type(paw_ij_type),target,intent(inout) :: paw_ij(my_natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
 type(pawrhoij_type),intent(in),optional :: electronpositron_pawrhoij(:)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
!Possible algos for PAW+U: 1=using occupation matrix n_i,,2=using PAW matrix rho_ij
 integer, parameter :: PAWU_ALGO_1=1,PAWU_ALGO_2=2
 integer, parameter :: PAWU_FLL=1,PAWU_AMF=2
 integer :: cplex_dij,iatom,iatom_tot,idij,ipositron,itypat,klmn,klmn1,lm_size,lmn2_size
 integer :: lpawu,my_comm_atom,my_comm_grid,natvshift_,ndij,nsploop,nsppol
 integer :: pawu_algo,pawu_dblec,qphase,usekden,usepawu,usexcnhat
 logical :: dij_available,dij_need,dij_prereq
 logical :: dij0_available,dij0_need,dij0_prereq
 logical :: dijexxc_available,dijexxc_need,dijexxc_prereq
 logical :: dijfock_available,dijfock_need,dijfock_prereq
 logical :: dijhartree_available,dijhartree_need,dijhartree_prereq
 logical :: dijhat_available,dijhat_need,dijhat_prereq
 logical :: dijhatfr_available,dijhatfr_need,dijhatfr_prereq
 logical :: dijnd_available,dijnd_need,dijnd_prereq
 logical :: dijso_available,dijso_need,dijso_prereq
 logical :: dijxc_available,dijxc_need,dijxc_prereq
 logical :: dijxchat_available,dijxchat_need,dijxchat_prereq
 logical :: dijxcval_available,dijxcval_need,dijxcval_prereq
 logical :: dijU_available,dijU_need,dijU_prereq
 logical :: has_nucdipmom,my_atmtab_allocated
 logical :: need_to_print,paral_atom,v_dijhat_allocated
 real(dp) :: hyb_mixing_,hyb_mixing_sr_
 character(len=500) :: msg
!arrays
 integer,pointer :: my_atmtab(:)
 logical,allocatable :: lmselect(:)
 real(dp),allocatable :: dij0(:),dijhartree(:)
 real(dp),allocatable :: dijhat(:,:),dijexxc(:,:),dijfock_cv(:,:),dijfock_vv(:,:),dijpawu(:,:)
 real(dp),allocatable :: dijnd(:,:),dijso(:,:),dijxc(:,:),dij_ep(:),dijxchat(:,:),dijxcval(:,:)
 real(dp),pointer :: v_dijhat(:,:),vpawu(:,:,:,:),vpawx(:,:,:)

! *************************************************************************

!------------------------------------------------------------------------
!----- Check consistency of arguments
!------------------------------------------------------------------------

!  === Check optional arguments ===

 hyb_mixing_   =zero ; if(present(hyb_mixing))    hyb_mixing_   =hyb_mixing
 hyb_mixing_sr_=zero ; if(present(hyb_mixing_sr)) hyb_mixing_sr_=hyb_mixing_sr

 natvshift_=0;if (present(natvshift)) natvshift_=natvshift
 if (natvshift_>0) then
   if ((.not.present(atvshift)).or.(.not.present(fatvshift))) then
     msg='when natvshift>0, atvshift and fatvshift arguments must be present!'
     MSG_BUG(msg)
   end if
 end if

 ipositron=0;if (present(electronpositron_calctype)) ipositron=electronpositron_calctype
 if (ipositron/=0) then
   if ((.not.present(electronpositron_pawrhoij)).or.&
&      (.not.present(electronpositron_lmselect))) then
     msg='ep_pawrhoij and ep_lmselect must be present for electron-positron calculations!'
     MSG_BUG(msg)
   end if
 end if

 has_nucdipmom=present(nucdipmom)

!  === Check complex character of arguments ===

 if (nspden==4.and.cplex==2) then
   msg='nspden=4 probably not compatible with cplex=2!'
   MSG_BUG(msg)
 end if
 if (my_natom>0) then
   if (paw_ij(1)%ndij==4.and.paw_ij(1)%cplex_dij/=2) then
     msg='invalid cplex size for Dij (4 Dij components)!'
     MSG_BUG(msg)
   end if
   if (paw_ij(1)%qphase/=paw_an(1)%cplex) then
     msg='paw_ij()%qphase and paw_an()%cplex must be equal!'
     MSG_BUG(msg)
   end if
   if (ipert<=0.and.paw_ij(1)%qphase/=1) then
     msg='qphase must be 1 for GS calculations!'
     MSG_BUG(msg)
   end if
   if (ipert>0.and.paw_ij(1)%qphase/=cplex) then
     msg='paw_ij()%qphase must be equal to cplex!'
     MSG_BUG(msg)
   end if
   if (paw_an(1)%has_vxcval>0.and.paw_an(1)%has_vxctau==2) then
     msg='kinetic energy density not available for vxc_val!'
     MSG_BUG(msg)
   end if
 end if

!------------------------------------------------------------------------
!----- Initializations
!------------------------------------------------------------------------

!Nothing to do for some perturbations (RF case)
 if (ipert==natom+1.or.ipert==natom+10) then
   do iatom=1,my_natom
     if (paw_ij(iatom)%has_dij==1) paw_ij(iatom)%dij=zero
     if (paw_ij(iatom)%has_dij0==1) paw_ij(iatom)%dij0=zero
     if (paw_ij(iatom)%has_dijfock==1) paw_ij(iatom)%dijfock=zero
     if (paw_ij(iatom)%has_dijhartree==1) paw_ij(iatom)%dijhartree=zero
     if (paw_ij(iatom)%has_dijxc==1) paw_ij(iatom)%dijxc=zero
     if (paw_ij(iatom)%has_dijhat==1) paw_ij(iatom)%dijhat=zero
     if (paw_ij(iatom)%has_dijso==1) paw_ij(iatom)%dijso=zero
     if (paw_ij(iatom)%has_dijU==1) paw_ij(iatom)%dijU=zero
     if (paw_ij(iatom)%has_dijexxc==1) paw_ij(iatom)%dijexxc=zero
     if (paw_ij(iatom)%has_dijxc_hat==1) paw_ij(iatom)%dijxc_hat=zero
     if (paw_ij(iatom)%has_dijxc_val==1) paw_ij(iatom)%dijxc_val=zero
   end do
   return
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!----- Various initializations
 nsppol=1;nsploop=1
 if (my_natom>0) then
   nsppol=paw_ij(1)%nsppol
   nsploop=nsppol;if (paw_ij(1)%ndij==4) nsploop=4
 end if
 usexcnhat=maxval(pawtab(1:ntypat)%usexcnhat)
 my_comm_grid=xmpi_comm_self;if (present(mpi_comm_grid)) my_comm_grid=mpi_comm_grid

!------ Select potential for Dij^hat computation
 v_dijhat_allocated=.false.
 if (my_natom>0) then
   if ((paw_ij(1)%has_dij==1).or.(paw_ij(1)%has_dijhat==1).or. &
&      (paw_ij(1)%has_dijhat==0.and.pawprtvol/=0)) then
     if (usexcnhat==0) then
       if (size(vxc,1)/=cplex*nfft.or.size(vxc,2)/=nspden) then
         msg='invalid size for vxc!'
         MSG_BUG(msg)
       end if
       LIBPAW_POINTER_ALLOCATE(v_dijhat,(cplex*nfft,nspden))
       v_dijhat_allocated=.true.
       !v_dijhat=vtrial-vxc
       do idij=1,nspden
         do klmn=1,cplex*nfft
           v_dijhat(klmn,idij)=vtrial(klmn,idij)-vxc(klmn,idij)
         end do
       end do
     else
       v_dijhat => vtrial
     end if
   end if
 end if

!------------------------------------------------------------------------
!----- Loop over atoms
!------------------------------------------------------------------------

 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

!  === Atom-dependent data ===

   itypat=paw_ij(iatom)%itypat
   cplex_dij=paw_ij(iatom)%cplex_dij
   qphase=paw_ij(iatom)%qphase
   lm_size=paw_an(iatom)%lm_size
   lmn2_size=paw_ij(iatom)%lmn2_size
   ndij=paw_ij(iatom)%ndij
   usepawu=pawtab(itypat)%usepawu
   pawu_algo=merge(PAWU_ALGO_1,PAWU_ALGO_2,ipert<=0.and.usepawu>=0)
   pawu_dblec=merge(PAWU_FLL,PAWU_AMF,abs(usepawu)==1.or.abs(usepawu)==4)
   usekden=merge(0,1,paw_an(iatom)%has_vxctau/=2)
   need_to_print=((abs(pawprtvol)>=1).and. &
&   (iatom_tot==1.or.iatom_tot==natom.or.pawprtvol<0))

!  === Determine which conditions and prerequisites are fulfilled for Dij ===

 if (my_natom>0) then
!  Total Dij: no condition ; no prerequisites
   dij_available=.true.;dij_prereq=.true.
!  Dij0: not available for RF ; need kij for the positron
   dij0_available=(ipert<=0);dij0_prereq=(ipositron/=1.or.pawtab(itypat)%has_kij==2)
!  DijFock:not available for RF, positron; only for Fock exact exch. ; Vxc_ex needed
   dijfock_available=(paw_ij(iatom)%has_dijfock>0.and.ipert<=0.and.ipositron/=1)
   dijfock_prereq=(paw_ij(iatom)%has_dijfock==2)
!  DijHartree: no condition ; no prerequisites
   dijhartree_available=.true.;dijhartree_prereq=.true.
!  DijXC: no condition ; Vxc needed
   dijxc_available=.true.
   dijxc_prereq=(paw_ij(iatom)%has_dijxc==2.or.paw_an(iatom)%has_vxc>0)
!  Dij^hat: no condition ; no prerequisites
   dijhat_available=.true.;dijhat_prereq=.true.
!  Dij^hat_FR: only for RF and when it was previously computed
   dijhatfr_available=(ipert>0.and.paw_ij(iatom)%has_dijfr==2) ; dijhatfr_prereq=.true.
!  DijND: not available for RF, requires non-zero nucdipmom
   dijnd_available=.false. ; dijnd_prereq=(cplex_dij==2)
   if (has_nucdipmom) dijnd_available=(ipert<=0.and.any(abs(nucdipmom(:,iatom))>tol8))
!  DijSO: not available for RF, positron; only for spin-orbit ; VHartree and Vxc needed
   dijso_available=(pawspnorb>0.and.ipert<=0.and.ipositron/=1)
   dijso_prereq=(paw_ij(iatom)%has_dijso==2.or.&
&               (paw_an(iatom)%has_vhartree>0.and.paw_an(iatom)%has_vxc>0))
!  DijU: not available for positron; only for LDA+U
   dijU_available=(pawtab(itypat)%usepawu/=0.and.ipositron/=1)
   dijU_prereq=(paw_ij(iatom)%has_dijU==2.or.paw_ij(iatom)%has_pawu_occ>0.or. &
&              (paw_ij(iatom)%has_dijU>0))
!  DijExxc: not available for RF, positron; only for local exact exch. ; Vxc_ex needed
   dijexxc_available=(pawtab(itypat)%useexexch/=0.and.ipert<=0.and.ipositron/=1)
   dijexxc_prereq=(paw_ij(iatom)%has_dijexxc==2.or.paw_ij(iatom)%has_exexch_pot>0)
!  DijXC^hat: not available for RF ; Vxc needed
   dijxchat_available=(ipert<=0)
   dijxchat_prereq=(paw_ij(iatom)%has_dijxc_hat==2.or.paw_an(iatom)%has_vxc>0)
!  DijXC_val: not available for RF ; Vxc_val needed
   dijxcval_available=(ipert<=0)
   dijxcval_prereq=(paw_ij(iatom)%has_dijxc_val==2.or.paw_an(iatom)%has_vxcval>0)
 end if

!  === Determine which parts of Dij have to be computed ===

   dij_need=.false.;dij0_need=.false.;dijexxc_need=.false.;dijfock_need=.false.
   dijhartree_need=.false.;dijhat_need=.false.;dijhatfr_need=.false.;
   dijso_need=.false.;dijU_need=.false.;dijxc_need=.false.;dijxchat_need=.false.
   dijxcval_need=.false.; dijnd_need=.false.

   if (dij_available) then
     if (paw_ij(iatom)%has_dij==1) then
       dij_need=.true.;paw_ij(iatom)%dij(:,:)=zero
     else if (paw_ij(iatom)%has_dij==0.and.need_to_print) then
       LIBPAW_ALLOCATE(paw_ij(iatom)%dij,(cplex_dij*qphase*lmn2_size,ndij))
       dij_need=.true.;paw_ij(iatom)%dij(:,:)=zero
       paw_ij(iatom)%has_dij=-1
     end if
   else if (paw_ij(iatom)%has_dij==1) then
     paw_ij(iatom)%dij=zero
   end if

   if (dij0_available) then
     if (paw_ij(iatom)%has_dij0==1) then
       dij0_need=.true.;paw_ij(iatom)%dij0(:)=zero
     else if (paw_ij(iatom)%has_dij0==0.and.need_to_print) then
       LIBPAW_ALLOCATE(paw_ij(iatom)%dij0,(lmn2_size))
       dij0_need=.true.;paw_ij(iatom)%dij0(:)=zero
       paw_ij(iatom)%has_dij0=-1
     end if
   else if (paw_ij(iatom)%has_dij0==1) then
     paw_ij(iatom)%dij0=zero
   end if

   if (dijfock_available) then
     if (paw_ij(iatom)%has_dijfock==1) then
       dijfock_need=.true.;paw_ij(iatom)%dijfock(:,:)=zero
     else if (paw_ij(iatom)%has_dijfock==0.and.need_to_print) then
       LIBPAW_ALLOCATE(paw_ij(iatom)%dijfock,(cplex_dij*lmn2_size,ndij))
       dijfock_need=.true.;paw_ij(iatom)%dijfock(:,:)=zero
       paw_ij(iatom)%has_dijfock=-1
     end if
   else if (paw_ij(iatom)%has_dijfock==1) then
     paw_ij(iatom)%dijfock=zero
   end if

   if (dijhartree_available) then
     if (paw_ij(iatom)%has_dijhartree==1) then
       dijhartree_need=.true.;paw_ij(iatom)%dijhartree(:)=zero
     else if (paw_ij(iatom)%has_dijhartree==0) then
       LIBPAW_ALLOCATE(paw_ij(iatom)%dijhartree,(qphase*lmn2_size))
       dijhartree_need=.true.;paw_ij(iatom)%dijhartree(:)=zero
       paw_ij(iatom)%has_dijhartree=-1
     end if
   else if (paw_ij(iatom)%has_dijhartree==1) then
     paw_ij(iatom)%dijhartree=zero
   end if

   if (dijxc_available) then
     if (paw_ij(iatom)%has_dijxc==1) then
       dijxc_need=.true.;paw_ij(iatom)%dijxc(:,:)=zero
     else if (paw_ij(iatom)%has_dijxc==0.and.need_to_print) then
       LIBPAW_ALLOCATE(paw_ij(iatom)%dijxc,(cplex_dij*qphase*lmn2_size,ndij))
       dijxc_need=.true.;paw_ij(iatom)%dijxc(:,:)=zero
       paw_ij(iatom)%has_dijxc=-1
     end if
   else if (paw_ij(iatom)%has_dijxc==1) then
     paw_ij(iatom)%dijxc=zero
   end if

   if (dijhat_available) then
     if (paw_ij(iatom)%has_dijhat==1) then
       dijhat_need=.true.;paw_ij(iatom)%dijhat(:,:)=zero
     else if (paw_ij(iatom)%has_dijhat==0.and.need_to_print) then
       LIBPAW_ALLOCATE(paw_ij(iatom)%dijhat,(cplex_dij*qphase*lmn2_size,ndij))
       dijhat_need=.true.;paw_ij(iatom)%dijhat(:,:)=zero
      paw_ij(iatom)%has_dijhat=-1
     end if
   else if (paw_ij(iatom)%has_dijhat==1) then
     paw_ij(iatom)%dijhat=zero
   end if

   if (dijnd_available) then
     if (paw_ij(iatom)%has_dijnd==1) then
       dijnd_need=.true.;paw_ij(iatom)%dijnd(:,:)=zero
     else if (paw_ij(iatom)%has_dijnd==0.and.need_to_print) then
       LIBPAW_ALLOCATE(paw_ij(iatom)%dijnd,(cplex_dij*lmn2_size,ndij))
       dijnd_need=.true.;paw_ij(iatom)%dijnd(:,:)=zero
       paw_ij(iatom)%has_dijnd=-1
     end if
   else if (paw_ij(iatom)%has_dijnd==1) then
     paw_ij(iatom)%dijnd=zero
   end if

   if (dijso_available) then
     if (paw_ij(iatom)%has_dijso==1) then
       dijso_need=.true.;paw_ij(iatom)%dijso(:,:)=zero
     else if (paw_ij(iatom)%has_dijso==0.and.need_to_print) then
       LIBPAW_ALLOCATE(paw_ij(iatom)%dijso,(cplex_dij*qphase*lmn2_size,ndij))
       dijso_need=.true.;paw_ij(iatom)%dijso(:,:)=zero
       paw_ij(iatom)%has_dijso=-1
     end if
   else if (paw_ij(iatom)%has_dijso==1) then
     paw_ij(iatom)%dijso=zero
   end if

   if (dijU_available) then
     if (paw_ij(iatom)%has_dijU==1) then
       dijU_need=.true.;paw_ij(iatom)%dijU(:,:)=zero
     else if (paw_ij(iatom)%has_dijU==0.and.need_to_print) then
       LIBPAW_ALLOCATE(paw_ij(iatom)%dijU,(cplex_dij*qphase*lmn2_size,ndij))
       dijU_need=.true.;paw_ij(iatom)%dijU(:,:)=zero
       paw_ij(iatom)%has_dijU=-1
     end if
   else if (paw_ij(iatom)%has_dijU==1) then
     paw_ij(iatom)%dijU=zero
   end if

   if (dijexxc_available.and.paw_ij(iatom)%has_dijexxc/=2) then
     if (paw_ij(iatom)%has_dijexxc==1) then
       dijexxc_need=.true.;paw_ij(iatom)%dijexxc(:,:)=zero
     else if (paw_ij(iatom)%has_dijexxc==0.and.need_to_print) then
       LIBPAW_ALLOCATE(paw_ij(iatom)%dijexxc,(cplex_dij*lmn2_size,ndij))
       dijexxc_need=.true.;paw_ij(iatom)%dijexxc(:,:)=zero
       paw_ij(iatom)%has_dijexxc=-1
     end if
   else if (paw_ij(iatom)%has_dijexxc==1) then
     paw_ij(iatom)%dijexxc=zero
   end if

   if (dijxchat_available) then
     if (paw_ij(iatom)%has_dijxc_hat==1) then
       dijxchat_need=.true.;paw_ij(iatom)%dijxc_hat(:,:)=zero
!      else if (paw_ij(iatom)%has_dijxc_hat==0.and.need_to_print) then
!      LIBPAW_ALLOCATE(paw_ij(iatom)%dijxc_hat,(cplex_dij*qphase*lmn2_size,ndij))
!      dijxchat_need=.true.;paw_ij(iatom)%dijxc_hat(:,:)=zero
!      paw_ij(iatom)%has_dijxc_hat=-1
     end if
   else if (paw_ij(iatom)%has_dijxc_hat==1) then
     paw_ij(iatom)%dijxc_hat=zero
   end if

   if (dijxcval_available) then
     if (paw_ij(iatom)%has_dijxc_val==1) then
       dijxcval_need=.true.;paw_ij(iatom)%dijxc_val(:,:)=zero
!      else if (paw_ij(iatom)%has_dijxc_val==0.and.need_to_print) then
!      LIBPAW_ALLOCATE(paw_ij(iatom)%dijxc_val,(cplex_dij*qphase*lmn2_size,ndij))
!      dijxcval_need=.true.;paw_ij(iatom)%dijxc_val(:,:)=zero
!      paw_ij(iatom)%has_dijxc_val=-1
     end if
   else if (paw_ij(iatom)%has_dijxc_val==1) then
     paw_ij(iatom)%dijxc_val=zero
   end if

!  === Print error messages if prerequisites are not fulfilled ===

   if (dij_need.and.(.not.dij_prereq)) then
     msg='Dij prerequisites missing!'
     MSG_BUG(msg)
   end if
   if (dij0_need.and.(.not.dij0_prereq)) then
     msg='Dij0 prerequisites missing!'
     MSG_BUG(msg)
   end if
   if (dijfock_need.and.(.not.dijfock_prereq)) then
     msg='DijFock prerequisites missing!'
     MSG_BUG(msg)
   end if

   if (dijhartree_need.and.(.not.dijhartree_prereq)) then
     msg='DijHartree prerequisites missing!'
     MSG_BUG(msg)
   end if
   if (dijxc_need.and.(.not.dijxc_prereq)) then
     msg='Dij^XC prerequisites missing!'
     MSG_BUG(msg)
   end if
   if (dijhat_need.and.(.not.dijhat_prereq)) then
     msg='Dij^hat prerequisites missing!'
     MSG_BUG(msg)
   end if
   if (dijhatfr_need.and.(.not.dijhatfr_prereq)) then
     msg='DijFR^hat prerequisites missing!'
     MSG_BUG(msg)
   end if
   if (dijnd_need.and.(.not.dijnd_prereq)) then
     msg='DijND prerequisites missing!'
     MSG_BUG(msg)
   end if
   if (dijso_need.and.(.not.dijso_prereq)) then
     msg='DijSO prerequisites missing!'
     MSG_BUG(msg)
   end if
   if (dijU_need.and.(.not.dijU_prereq)) then
     msg='DijU prerequisites missing!'
     MSG_BUG(msg)
   end if
   if (dijexxc_need.and.(.not.dijexxc_prereq)) then
     msg='DijExcc prerequisites missing!'
     MSG_BUG(msg)
   end if
   if (dijxchat_need.and.(.not.dijxchat_prereq)) then
     msg='DijXC^hat prerequisites missing!'
     MSG_BUG(msg)
   end if
   if (dijxcval_need.and.(.not.dijxcval_prereq)) then
     msg='DijXC_val prerequisites missing!'
     MSG_BUG(msg)
   end if

!  ------------------------------------------------------------------------
!  ----------- Add atomic Dij0 to Dij
!  ------------------------------------------------------------------------

   if ((dij0_need.or.dij_need).and.dij0_available) then

     LIBPAW_ALLOCATE(dij0,(lmn2_size))
!    ===== Dij0 already computed
     if (paw_ij(iatom)%has_dij0==2) then
       dij0(:)=paw_ij(iatom)%dij0(:)
     else
!    ===== Need to compute Dij0
       dij0(:)=pawtab(itypat)%dij0(:)
       if (ipositron==1) dij0(:)=two*pawtab(itypat)%kij(:)-dij0(:)
       if (pawu_algo==PAWU_ALGO_2.and.pawu_dblec==PAWU_FLL) dij0(:)=dij0(:)+pawtab(itypat)%euij_fll(:)
       if (dij0_need) paw_ij(iatom)%dij0(:)=dij0(:)
     end if

     if (dij_need) then
       do idij=1,min(nsploop,2)
         klmn1=1
         do klmn=1,lmn2_size
           paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+dij0(klmn)
           klmn1=klmn1+cplex_dij
         end do
       end do
     end if
     LIBPAW_DEALLOCATE(dij0)
   end if

!  ------------------------------------------------------------------------
!  ------------------------------------------------------------------------
!  ----------- Add Dij_{Fock exact-exchange} to Dij
!  ------------------------------------------------------------------------

   if ((dijfock_need.or.dij_need).and.dijfock_available) then

!    ===== DijFock already computed
     if (paw_ij(iatom)%has_dijfock==2) then
       if (dij_need) paw_ij(iatom)%dij(1:cplex_dij*lmn2_size,:)= &
&                    paw_ij(iatom)%dij(1:cplex_dij*lmn2_size,:) &
&                   +paw_ij(iatom)%dijfock(1:cplex_dij*lmn2_size,:)

     else

!    ===== Need to compute DijFock
       LIBPAW_ALLOCATE(dijfock_vv,(cplex_dij*lmn2_size,ndij))
       LIBPAW_ALLOCATE(dijfock_cv,(cplex_dij*lmn2_size,ndij))
       dijfock_vv(:,:)=zero ; dijfock_cv(:,:)=zero
!      Exact exchange is evaluated for electrons only
       if (ipositron/=1) then
         call pawdijfock(dijfock_vv,dijfock_cv,cplex_dij,qphase,hyb_mixing_,hyb_mixing_sr_, &
&                        ndij,pawrhoij(iatom),pawtab(itypat))
       end if
       if (dijfock_need) paw_ij(iatom)%dijfock(:,:)=dijfock_vv(:,:)+dijfock_cv(:,:)
       if (dij_need) paw_ij(iatom)%dij(1:cplex_dij*lmn2_size,:)= &
&                    paw_ij(iatom)%dij(1:cplex_dij*lmn2_size,:) &
&                   +dijfock_vv(1:cplex_dij*lmn2_size,:)+dijfock_cv(1:cplex_dij*lmn2_size,:)
       LIBPAW_DEALLOCATE(dijfock_vv)
       LIBPAW_DEALLOCATE(dijfock_cv)
     end if
   end if

!  ----------- Add Dij_Hartree to Dij
!  ------------------------------------------------------------------------

   if ((dijhartree_need.or.dij_need).and.dijhartree_available) then

     LIBPAW_ALLOCATE(dijhartree,(qphase*lmn2_size))
!    ===== DijHartree already computed
     if (paw_ij(iatom)%has_dijhartree==2) then
       dijhartree(:)=paw_ij(iatom)%dijhartree(:)
     else
!    ===== Need to compute DijHartree
       if (ipositron/=1) then
         call pawdijhartree(dijhartree,qphase,nspden,pawrhoij(iatom),pawtab(itypat))
       else
         dijhartree(:)=zero
       end if
       if (ipositron/=0) then
         LIBPAW_ALLOCATE(dij_ep,(qphase*lmn2_size))
         call pawdijhartree(dij_ep,qphase,nspden,electronpositron_pawrhoij(iatom),pawtab(itypat))
         dijhartree(:)=dijhartree(:)-dij_ep(:)
         LIBPAW_DEALLOCATE(dij_ep)
       end if
       if (dijhartree_need) paw_ij(iatom)%dijhartree(:)=dijhartree(:)
     end if

     if (dij_need) then
       do idij=1,min(nsploop,2)
         klmn1=1
         do klmn=1,qphase*lmn2_size
           paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+dijhartree(klmn)
           klmn1=klmn1+cplex_dij
         end do
       end do
     end if

     LIBPAW_DEALLOCATE(dijhartree)
   end if

!  ------------------------------------------------------------------------
!  ----------- Add Dij_xc to Dij
!  ------------------------------------------------------------------------

   if ((dijxc_need.or.dij_need).and.dijxc_available) then

!    ===== Dijxc already computed
     if (paw_ij(iatom)%has_dijxc==2) then
       if (dij_need) paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)+paw_ij(iatom)%dijxc(:,:)
     else

!    ===== Need to compute DijXC
       LIBPAW_ALLOCATE(dijxc,(cplex_dij*qphase*lmn2_size,ndij))
       if (pawxcdev/=0) then
         LIBPAW_ALLOCATE(lmselect,(lm_size))
         lmselect(:)=paw_an(iatom)%lmselect(:)
         if (ipositron/=0) lmselect(:)=(lmselect(:).or.electronpositron_lmselect(1:lm_size,iatom))
         call pawdijxcm(dijxc,cplex_dij,qphase,lmselect,ndij,nspden,nsppol,pawang,&
&                       pawrad(itypat),pawtab(itypat),paw_an(iatom)%vxc1,&
&                       paw_an(iatom)%vxct1,usexcnhat)
         LIBPAW_DEALLOCATE(lmselect)
       else
         call pawdijxc(dijxc,cplex_dij,qphase,ndij,nspden,nsppol,&
&                      pawang,pawrad(itypat),pawtab(itypat),paw_an(iatom)%vxc1,&
&                      paw_an(iatom)%vxct1,usexcnhat,usekden,&
&                      vxctau1=paw_an(iatom)%vxctau1,vxcttau1=paw_an(iatom)%vxcttau1)
       end if
       if (dijxc_need) paw_ij(iatom)%dijxc(:,:)=dijxc(:,:)
       if (dij_need) paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)+dijxc(:,:)
       LIBPAW_DEALLOCATE(dijxc)
     end if

   end if

!  ------------------------------------------------------------------------
!  ----------- Add Dij_hat to Dij
!  ------------------------------------------------------------------------

   if ((dijhat_need.or.dij_need).and.dijhat_available) then

!    ===== Dijhat already computed
     if (paw_ij(iatom)%has_dijhat==2) then
       if (dij_need) paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)+paw_ij(iatom)%dijhat(:,:)
     else

!    ===== Need to compute Dijhat
       LIBPAW_ALLOCATE(dijhat,(cplex_dij*qphase*lmn2_size,ndij))
       call pawdijhat(dijhat,cplex_dij,qphase,gprimd,iatom_tot,&
&                     natom,ndij,nfft,nfftot,nspden,nsppol,pawang,pawfgrtab(iatom),&
&                     pawtab(itypat),v_dijhat,qphon,ucvol,xred,mpi_comm_grid=my_comm_grid)
       if (dijhat_need) paw_ij(iatom)%dijhat(:,:)=dijhat(:,:)
       if (dij_need) paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)+dijhat(:,:)
       LIBPAW_DEALLOCATE(dijhat)
     end if

!    ===== RF: add frozen part of 1st-order Dij
     if (dijhatfr_available) then
       do idij=1,nsploop
         if (dij_need) paw_ij(iatom)%dij(:,idij)=paw_ij(iatom)%dij(:,idij) &
&                                               +paw_ij(iatom)%dijfr(:,idij)
         if (dijhat_need) paw_ij(iatom)%dijhat(:,idij)=paw_ij(iatom)%dijhat(:,idij) &
&                                                     +paw_ij(iatom)%dijfr(:,idij)
       end do
     end if

   end if

!  ------------------------------------------------------------------------
!  ----------- Add Dij nuclear dipole moments to Dij
!  ------------------------------------------------------------------------

   if ((dijnd_need.or.dij_need).and.dijnd_available) then

!    ===== Dijnd already computed
     if (paw_ij(iatom)%has_dijnd==2) then
       if (dij_need) paw_ij(iatom)%dij(1:cplex_dij*lmn2_size,:)= &
&                    paw_ij(iatom)%dij(1:cplex_dij*lmn2_size,:) &
&                   +paw_ij(iatom)%dijnd(1:cplex_dij*lmn2_size,:)
     else

!    ===== Need to compute Dijnd
       LIBPAW_ALLOCATE(dijnd,(cplex_dij*lmn2_size,ndij))
       call pawdijnd(dijnd,cplex_dij,ndij,nucdipmom(:,iatom),pawrad(itypat),pawtab(itypat))
       if (dijnd_need) paw_ij(iatom)%dijnd(:,:)=dijnd(:,:)
       if (dij_need) paw_ij(iatom)%dij(1:cplex_dij*lmn2_size,:)= &
&                    paw_ij(iatom)%dij(1:cplex_dij*lmn2_size,:) &
&                   +dijnd(1:cplex_dij*lmn2_size,:)
       LIBPAW_DEALLOCATE(dijnd)
     end if

   end if


!  ------------------------------------------------------------------------
!  ----------- Add Dij spin-orbit to Dij
!  ------------------------------------------------------------------------

   if ((dijso_need.or.dij_need).and.dijso_available) then

!    ===== DijSO already computed
     if (paw_ij(iatom)%has_dijso==2) then
       if (dij_need) paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)+paw_ij(iatom)%dijso(:,:)
     else

!    ===== Need to compute DijSO
       LIBPAW_ALLOCATE(dijso,(cplex_dij*qphase*lmn2_size,ndij))
       call pawdijso(dijso,cplex_dij,qphase,ndij,nspden,&
&                    pawang,pawrad(itypat),pawtab(itypat),pawxcdev,spnorbscl,&
&                    paw_an(iatom)%vh1,paw_an(iatom)%vxc1)
       if (dijso_need) paw_ij(iatom)%dijso(:,:)=dijso(:,:)
       if (dij_need) paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)+dijso(:,:)
       LIBPAW_DEALLOCATE(dijso)
     end if

   end if

!  ------------------------------------------------------------------------
!  ----------- Add Dij_{LDA+U} to Dij
!  ------------------------------------------------------------------------

   if ((dijU_need.or.dij_need).and.dijU_available) then

!    ===== DijU already computed
     if (paw_ij(iatom)%has_dijU==2) then
       if (dij_need) paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)+paw_ij(iatom)%dijU(:,:)
     else

!    ===== Need to compute DijU
       LIBPAW_ALLOCATE(dijpawu,(cplex_dij*qphase*lmn2_size,ndij))
       if (pawu_algo==PAWU_ALGO_2) then
         call pawdiju_euijkl(dijpawu,cplex_dij,qphase,ndij,pawrhoij(iatom),pawtab(itypat))
       else
         lpawu=pawtab(itypat)%lpawu
         LIBPAW_POINTER_ALLOCATE(vpawu,(cplex_dij,lpawu*2+1,lpawu*2+1,ndij))
         if (usepawu>=10) vpawu=zero ! if dmft, do not apply U in LDA+U
         if (usepawu< 10) then
           call pawpupot(cplex_dij,ndij,paw_ij(iatom)%noccmmp,paw_ij(iatom)%nocctot,&
&                        pawprtvol,pawtab(itypat),vpawu)
         end if
         if (natvshift_==0) then
           call pawdiju(dijpawu,cplex_dij,qphase,ndij,nsppol,pawtab(itypat),vpawu)
         else
           call pawdiju(dijpawu,cplex_dij,qphase,ndij,nsppol,pawtab(itypat),vpawu,&
&                       natvshift=natvshift_,atvshift=atvshift(:,:,iatom_tot),&
&                       fatvshift=fatvshift)
         end if
         LIBPAW_POINTER_DEALLOCATE(vpawu)
       end if
       if (dijU_need) paw_ij(iatom)%dijU(:,:)=dijpawu(:,:)
       if (dij_need) paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)+dijpawu(:,:)
       LIBPAW_DEALLOCATE(dijpawu)
     end if

   end if

!  ------------------------------------------------------------------------
!  ----------- Add Dij_{local exact-exchange} to Dij
!  ------------------------------------------------------------------------

   if ((dijexxc_need.or.dij_need).and.dijexxc_available) then

!    ===== DijEXXC already computed
     if (paw_ij(iatom)%has_dijexxc==2) then
       if (dij_need) paw_ij(iatom)%dij(1:cplex_dij*lmn2_size,:)= &
&                    paw_ij(iatom)%dij(1:cplex_dij*lmn2_size,:) &
&                   +paw_ij(iatom)%dijexxc(1:cplex_dij*lmn2_size,:)
     else

!    ===== Need to compute DijEXXC
       LIBPAW_ALLOCATE(dijexxc,(cplex_dij*lmn2_size,ndij))
       if (pawxcdev/=0) then
         if (paw_ij(iatom)%has_exexch_pot/=2) then
           LIBPAW_POINTER_ALLOCATE(vpawx,(1,lmn2_size,ndij))
           call pawxpot(ndij,pawprtvol,pawrhoij(iatom),pawtab(itypat),vpawx)
         else
           vpawx=>paw_ij(iatom)%vpawx
         end if
         LIBPAW_ALLOCATE(lmselect,(lm_size))
         lmselect(:)=paw_an(iatom)%lmselect(:)
         if (ipositron/=0) lmselect(:)=(lmselect(:).or.electronpositron_lmselect(1:lm_size,iatom))
         call pawdijexxc(dijexxc,cplex_dij,qphase,lmselect,ndij,nspden,nsppol,&
&             pawang,pawrad(itypat),pawtab(itypat),vpawx,paw_an(iatom)%vxc_ex)
         LIBPAW_DEALLOCATE(lmselect)
         if (paw_ij(iatom)%has_exexch_pot/=2) then
            LIBPAW_POINTER_DEALLOCATE(vpawx)
         end if
         if (dijexxc_need) paw_ij(iatom)%dijexxc(:,:)=dijexxc(:,:)
         if (dij_need) paw_ij(iatom)%dij(1:cplex_dij*lmn2_size,:)= &
&                      paw_ij(iatom)%dij(1:cplex_dij*lmn2_size,:) &
&                     +dijexxc(1:cplex_dij*lmn2_size,:)
         LIBPAW_DEALLOCATE(dijexxc)
       end if
     end if

   end if

!  ------------------------------------------------------------------------
!  ----------- Add Dij background contribution to the total Dij
!  ------------------------------------------------------------------------

   if (dij_need.and.pawtab(itypat)%usepotzero==1 ) then
     do idij=1,min(nsploop,2)
       klmn1=1
       do klmn=1,lmn2_size
         paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+pawtab(itypat)%gammaij(klmn)*charge/ucvol
         klmn1=klmn1+cplex_dij*qphase
       end do
     end do
   end if


!  ------------------------------------------------------------------------
!  ----------- Compute Dijxc_hat
!  ------------------------------------------------------------------------

   if (dijxchat_need) then

     if (usexcnhat/=0) then
       LIBPAW_ALLOCATE(dijxchat,(cplex_dij*lmn2_size,ndij))
       call pawdijhat(dijxchat,cplex_dij,1,gprimd,iatom_tot,&
&                     natom,ndij,nfft,nfftot,nspden,nsppol,pawang,pawfgrtab(iatom),&
&                     pawtab(itypat),vxc,qphon,ucvol,xred,mpi_comm_grid=my_comm_grid)
       paw_ij(iatom)%dijxc_hat(1:cplex_dij*lmn2_size,:)=dijxchat(1:cplex_dij*lmn2_size,:)
       LIBPAW_DEALLOCATE(dijxchat)

     else ! usexcnhat=0
       paw_ij(iatom)%dijxc_hat=zero
     end if

   end if

!  ------------------------------------------------------------------------
!  ----------- Compute Dijxc_val
!  ------------------------------------------------------------------------

   if (dijxcval_need) then

     LIBPAW_ALLOCATE(dijxcval,(cplex_dij*lmn2_size,ndij))
!    Note that usexcnhat=0 for this call (no compensation term)
     if (pawxcdev/=0) then
       LIBPAW_ALLOCATE(lmselect,(lm_size))
       lmselect(:)=paw_an(iatom)%lmselect(:)
       if (ipositron/=0) lmselect(:)=(lmselect(:).or.electronpositron_lmselect(1:lm_size,iatom))
       call pawdijxcm(dijxcval,cplex_dij,1,lmselect,ndij,nspden,nsppol,&
&                     pawang,pawrad(itypat),pawtab(itypat),paw_an(iatom)%vxc1_val,&
&                     paw_an(iatom)%vxct1_val,0)
       LIBPAW_DEALLOCATE(lmselect)
     else
       call pawdijxc(dijxcval,cplex_dij,1,ndij,nspden,nsppol,&
&                    pawang,pawrad(itypat),pawtab(itypat),paw_an(iatom)%vxc1_val,&
&                    paw_an(iatom)%vxct1_val,0,0)
     end if
     paw_ij(iatom)%dijxc_val(1:cplex_dij*lmn2_size,:)=dijxcval(1:cplex_dij*lmn2_size,:)
     LIBPAW_DEALLOCATE(dijxcval)

   end if

!  ------------------------------------------------------------------------

!  Update some flags
   if (dij_need.and.paw_ij(iatom)%has_dij>=1) paw_ij(iatom)%has_dij=2
   if (dij0_need.and.paw_ij(iatom)%has_dij0>=1) paw_ij(iatom)%has_dij0=2
   if (dijfock_need.and.paw_ij(iatom)%has_dijfock>=1) paw_ij(iatom)%has_dijfock=2

   if (dijhartree_need.and.paw_ij(iatom)%has_dijhartree>=1) paw_ij(iatom)%has_dijhartree=2
   if (dijxc_need.and.paw_ij(iatom)%has_dijxc>=1) paw_ij(iatom)%has_dijxc=2
   if (dijhat_need.and.paw_ij(iatom)%has_dijhat>=1) paw_ij(iatom)%has_dijhat=2
   if (dijnd_need.and.paw_ij(iatom)%has_dijnd>=1) paw_ij(iatom)%has_dijnd=2
   if (dijso_need.and.paw_ij(iatom)%has_dijso>=1) paw_ij(iatom)%has_dijso=2
   if (dijU_need.and.paw_ij(iatom)%has_dijU>=1) paw_ij(iatom)%has_dijU=2
   if (dijexxc_need.and.paw_ij(iatom)%has_dijexxc>=1) paw_ij(iatom)%has_dijexxc=2
   if (dijxchat_need.and.paw_ij(iatom)%has_dijxc_hat>=1) paw_ij(iatom)%has_dijxc_hat=2
   if (dijxcval_need.and.paw_ij(iatom)%has_dijxc_val>=1) paw_ij(iatom)%has_dijxc_val=2

!End loop over atoms
 end do ! iatom

!------------------------------------------------------------------------

!Final printing
 if (paral_atom) then
   call paw_ij_print(paw_ij,unit=std_out,pawprtvol=pawprtvol,pawspnorb=pawspnorb,&
&   comm_atom=my_comm_atom,mpi_atmtab=my_atmtab,natom=natom,&
&   mode_paral='PERS',enunit=enunit,ipert=ipert)
 else
   call paw_ij_print(paw_ij,unit=std_out,pawprtvol=pawprtvol,pawspnorb=pawspnorb,&
&   mode_paral='COLL',enunit=enunit,ipert=ipert)
 end if

!Free temporary storage
 if (v_dijhat_allocated) then
   LIBPAW_POINTER_DEALLOCATE(v_dijhat)
 end if
 do iatom=1,my_natom
   if (paw_ij(iatom)%has_dij0==-1) then
     LIBPAW_DEALLOCATE(paw_ij(iatom)%dij0)
     paw_ij(iatom)%has_dij0=0
   end if
   if (paw_ij(iatom)%has_dijfock==-1) then
     LIBPAW_DEALLOCATE(paw_ij(iatom)%dijfock)
     paw_ij(iatom)%has_dijfock=0
   end if

   if (paw_ij(iatom)%has_dijhartree==-1) then
     LIBPAW_DEALLOCATE(paw_ij(iatom)%dijhartree)
     paw_ij(iatom)%has_dijhartree=0
   end if
   if (paw_ij(iatom)%has_dijxc==-1) then
     LIBPAW_DEALLOCATE(paw_ij(iatom)%dijxc)
     paw_ij(iatom)%has_dijxc=0
   end if
   if (paw_ij(iatom)%has_dijhat==-1) then
     LIBPAW_DEALLOCATE(paw_ij(iatom)%dijhat)
     paw_ij(iatom)%has_dijhat=0
   end if
   if (paw_ij(iatom)%has_dijfr==-1) then
     LIBPAW_DEALLOCATE(paw_ij(iatom)%dijfr)
     paw_ij(iatom)%has_dijfr=0
   end if
   if (paw_ij(iatom)%has_dijso==-1) then
     LIBPAW_DEALLOCATE(paw_ij(iatom)%dijso)
     paw_ij(iatom)%has_dijso=0
   end if
   if (paw_ij(iatom)%has_dijU==-1) then
     LIBPAW_DEALLOCATE(paw_ij(iatom)%dijU)
     paw_ij(iatom)%has_dijU=0
   end if
   if (paw_ij(iatom)%has_dijexxc==-1) then
     LIBPAW_DEALLOCATE(paw_ij(iatom)%dijexxc)
     paw_ij(iatom)%has_dijexxc=0
   end if
   if (paw_ij(iatom)%has_dijxc_hat==-1) then
     LIBPAW_DEALLOCATE(paw_ij(iatom)%dijxc_hat)
     paw_ij(iatom)%has_dijxc_hat=0
   end if
   if (paw_ij(iatom)%has_dijxc_val==-1) then
     LIBPAW_DEALLOCATE(paw_ij(iatom)%dijxc_val)
     paw_ij(iatom)%has_dijxc_val=0
   end if
 end do

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine pawdij
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawdijhartree
!! NAME
!! pawdijhartree
!!
!! FUNCTION
!! Compute the Hartree contribution to the PAW pseudopotential strength Dij
!! (for one atom only)
!!
!! INPUTS
!!  qphase=2 if dij contains a exp(-i.q.r) phase (as in the q<>0 RF case), 1 if not
!!  nspden=number of spin density components
!!  pawrhoij <type(pawrhoij_type)>= paw rhoij occupancies (and related data) for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data, for current atom
!!
!! OUTPUT
!!  dijhartree(qphase*lmn2_size)=  D_ij^Hartree terms
!!    When a exp(-i.q.r) phase is included (qphase=2):
!!      dij(1:lmn2_size,:)
!!          contains the real part of the phase, i.e. D_ij*cos(q.r)
!!      dij(lmn2_size+1:2*lmn2_size,:)
!!          contains the imaginary part of the phase, i.e. D_ij*sin(q.r)
!!
!! PARENTS
!!      m_pawdij,pawdenpot,pawdfptenergy
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijhartree(dijhartree,qphase,nspden,pawrhoij,pawtab)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: nspden,qphase
!arrays
 real(dp),intent(out) :: dijhartree(:)
 type(pawrhoij_type),intent(in) :: pawrhoij
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: cplex_rhoij,iq,iq0_dij,iq0_rhoij,irhoij,ispden,jrhoij,klmn,klmn1,lmn2_size,nspdiag
 real(dp) :: ro
 character(len=500) :: msg
!arrays

! *************************************************************************

!Check data consistency
 if (size(dijhartree,1)/=qphase*pawtab%lmn2_size) then
   msg='invalid size for DijHartree!'
   MSG_BUG(msg)
 end if
 if (pawrhoij%qphase<qphase) then
   msg='pawrhoij%qphase must be >=qphase!'
   MSG_BUG(msg)
 end if

!Initialization
 dijhartree=zero
 lmn2_size=pawrhoij%lmn2_size
 cplex_rhoij=pawrhoij%cplex_rhoij

!Loop over (diagonal) spin-components
 nspdiag=1;if (nspden==2) nspdiag=2
 do ispden=1,nspdiag

   !Loop over phase exp(iqr) phase real/imaginary part
   do iq=1,qphase
     !First loop: we store the real part in dij(1 -> lmn2_size)
     !2nd loop: we store the imaginary part in dij(lmn2_size+1 -> 2*lmn2_size)
     iq0_dij=merge(0,lmn2_size,iq==1)
     iq0_rhoij=cplex_rhoij*iq0_dij

     !Loop over rhoij elements
     jrhoij=iq0_rhoij+1
     do irhoij=1,pawrhoij%nrhoijsel
       klmn=pawrhoij%rhoijselect(irhoij)

       ro=pawrhoij%rhoijp(jrhoij,ispden)*pawtab%dltij(klmn)

       !Diagonal k=l
       dijhartree(iq0_dij+klmn)=dijhartree(iq0_dij+klmn)+ro*pawtab%eijkl(klmn,klmn)

       !k<=l
       do klmn1=1,klmn-1
         dijhartree(iq0_dij+klmn1)=dijhartree(iq0_dij+klmn1)+ro*pawtab%eijkl(klmn1,klmn)
       end do

       !k>l
       do klmn1=klmn+1,lmn2_size
         dijhartree(iq0_dij+klmn1)=dijhartree(iq0_dij+klmn1)+ro*pawtab%eijkl(klmn,klmn1)
       end do

       jrhoij=jrhoij+cplex_rhoij
     end do !End loop over rhoij

   end do !End loop over q phase

 end do !End loop over spin

end subroutine pawdijhartree
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawdijfock
!! NAME
!! pawdijfock
!!
!! FUNCTION
!! Compute Fock exact-exchange contribution(s) to the PAW pseudopotential strength Dij
!! (for one atom only)
!!
!! INPUTS
!!  hyb_mixing=hybrid mixing coefficient for the Fock contribution
!!  hyb_mixing_sr=hybrid mixing coefficient for the short-range Fock contribution
!!  qphase=2 if dij contains a exp(-i.q.r) phase (as in the q<>0 RF case), 1 if not
!!  ndij= number of spin components dor Fdij
!!  pawrhoij <type(pawrhoij_type)>= paw rhoij occupancies (and related data) for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data, for current atom
!!
!! OUTPUT
!!  dijfock_vv(qphase*lmn2_size,ndij)=  D_ij^fock terms for valence-valence interactions
!!  dijfock_cv(qphase*lmn2_size,ndij)=  D_ij^fock terms for core-valence interactions
!!    When a exp(-i.q.r) phase is included (qphase=2):
!!      dij(1:lmn2_size,:)
!!          contains the real part of the phase, i.e. D_ij*cos(q.r)
!!      dij(lmn2_size+1:2*lmn2_size,:)
!!          contains the imaginary part of the phase, i.e. D_ij*sin(q.r)
!!
!!  NOTES:
!!   WARNING: What follows has been tested only for cases where nsppol=1 and 2, nspden=1 and 2 with nspinor=1.
!!
!! PARENTS
!!      m_pawdij,pawdenpot
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijfock(dijfock_vv,dijfock_cv,cplex_dij,qphase,hyb_mixing,hyb_mixing_sr,ndij,pawrhoij,pawtab)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex_dij,ndij,qphase
 real(dp),intent(in) :: hyb_mixing,hyb_mixing_sr
!arrays
 real(dp),intent(out) :: dijfock_vv(:,:),dijfock_cv(:,:)
 type(pawrhoij_type),intent(in) :: pawrhoij
 type(pawtab_type),intent(in),target :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: cplex_rhoij,iq,iq0_dij,iq0_rhoij,ispden,irhokl,jrhokl,ilmn_i,jlmn_j,ilmn_k,jlmn_l
 integer :: klmn_kl,klmn_ij,klmn_il,klmn_kj,klmn1,nsp,lmn2_size
 real(dp) :: ro,dij_up,dij_dn,dij_updn_r,dij_updn_i
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: dijfock_vv_tmp(:,:)
 real(dp),pointer :: eijkl(:,:)

! *************************************************************************

!Useful data
 lmn2_size=pawtab%lmn2_size
 cplex_rhoij=pawrhoij%cplex_rhoij

!Check data consistency
 if (size(dijfock_vv,1)/=qphase*cplex_dij*lmn2_size.or.size(dijfock_vv,2)/=ndij) then
   msg='invalid sizes for Dijfock_vv!'
   MSG_BUG(msg)
 end if
 if (size(dijfock_cv,1)/=qphase*cplex_dij*lmn2_size.or.size(dijfock_cv,2)/=ndij) then
   msg='invalid sizes for Dijfock_cv!'
   MSG_BUG(msg)
 end if
 if (pawrhoij%qphase<qphase) then
   msg='pawrhoij%qphase must be >=qphase!'
   MSG_BUG(msg)
 end if
 if (ndij==4.and.cplex_dij==1) then
   msg='When ndij=4, cplex_dij must be =2!'
   MSG_BUG(msg)
 end if

 if (abs(hyb_mixing)>tol8 .and. abs(hyb_mixing_sr)>tol8) then
   msg='invalid hybrid functional'
   MSG_BUG(msg)
 else
   if (abs(hyb_mixing)>tol8) then
     eijkl => pawtab%eijkl
   else if (abs(hyb_mixing_sr)>tol8) then
     eijkl => pawtab%eijkl_sr
   end if
 end if

!Init memory
 dijfock_vv=zero ; dijfock_cv=zero

! ===== Valence-valence contribution =====

 nsp=pawrhoij%nsppol;if (pawrhoij%nspden==4) nsp=4
 LIBPAW_ALLOCATE(dijfock_vv_tmp,(lmn2_size,nsp))
 dijfock_vv_tmp=zero

!Loop over phase exp(iqr) phase real/imaginary part
 do iq=1,qphase
   !First loop: we store the real part in dij(1 -> lmn2_size)
   !2nd loop: we store the imaginary part in dij(lmn2_size+1 -> 2*lmn2_size)
   iq0_dij=merge(0,lmn2_size,iq==1) ; iq0_rhoij=cplex_rhoij*iq0_dij

   !Loop over spin components
   do ispden=1,nsp

     !Loop on the non-zero elements rho_kl
     jrhokl=iq0_rhoij+1
     do irhokl=1,pawrhoij%nrhoijsel
       klmn_kl=pawrhoij%rhoijselect(irhokl)
       ilmn_k=pawtab%indklmn(7,klmn_kl)
       jlmn_l=pawtab%indklmn(8,klmn_kl)

       ro=pawrhoij%rhoijp(jrhokl,ispden)*pawtab%dltij(klmn_kl)

       !Contribution to the element (k,l) of dijfock
       dijfock_vv_tmp(klmn_kl,ispden)=dijfock_vv_tmp(klmn_kl,ispden)-ro*eijkl(klmn_kl,klmn_kl)

       !Contribution to the element (i,j) of dijfock with (i,j) < (k,l)
       !  We remind that i<j and k<l by construction
       do klmn_ij=1,klmn_kl-1
         ilmn_i=pawtab%indklmn(7,klmn_ij)
         jlmn_j=pawtab%indklmn(8,klmn_ij)
         !In this case, i < l
         klmn_il=jlmn_l*(jlmn_l-1)/2+ilmn_i
         !For (k,j), we compute index of (k,j) or index of (j,k)
         if (ilmn_k>jlmn_j) then
           klmn_kj=ilmn_k*(ilmn_k-1)/2+jlmn_j
         else
           klmn_kj=jlmn_j*(jlmn_j-1)/2+ilmn_k
         end if
         dijfock_vv_tmp(klmn_ij,ispden)=dijfock_vv_tmp(klmn_ij,ispden)-ro*eijkl(klmn_il,klmn_kj)
       end do

       !Contribution to the element (i,j) of dijfock with (i,j) > (k,l)
       !  We remind that i<j and k<l by construction
       do klmn_ij=klmn_kl+1,lmn2_size
         ilmn_i=pawtab%indklmn(7,klmn_ij)
         jlmn_j=pawtab%indklmn(8,klmn_ij)
         !In this case, k < j
         klmn_kj=jlmn_j*(jlmn_j-1)/2+ilmn_k
         !For (i,l), we compute index of (i,l) or index of (l,i)
         if (ilmn_i>jlmn_l) then
           klmn_il=ilmn_i*(ilmn_i-1)/2+jlmn_l
         else
           klmn_il=jlmn_l*(jlmn_l-1)/2+ilmn_i
         end if
         dijfock_vv_tmp(klmn_ij,ispden)=dijfock_vv_tmp(klmn_ij,ispden)-ro*eijkl(klmn_kj,klmn_il)
       end do

       jrhokl=jrhokl+cplex_rhoij
     end do !End loop over rhoij

   end do !ispden

!  Regular case: copy spin component into Dij
   if (ndij/=4.or.nsp/=4) then
     do ispden=1,nsp
       klmn1=iq0_dij+1
       do klmn_ij=1,lmn2_size
         dijfock_vv(klmn1,ispden)=dijfock_vv_tmp(klmn_ij,ispden)
         klmn1=klmn1+cplex_dij
       end do
     end do
!    Antiferro case: copy up component into down one
     if (ndij==2.and.nsp==1) then
       klmn1=iq0_dij+1
       do klmn_ij=1,lmn2_size
         dijfock_vv(klmn1,2)=dijfock_vv_tmp(klmn_ij,1)
         klmn1=klmn1+cplex_dij
       end do
     end if
   else
   !Non-collinear: from (rhoij,m_ij) to rhoij^(alpha,beta)
   !rhoij=  (rhoij^11+rhoij^22)
   !mij_x=  (rhoij^12+rhoij^21)
   !mij_y=i.(rhoij^12+rhoij^21)
   !mij_z=  (rhoij^11-rhoij^22)
     klmn1=iq0_dij+1
     do klmn_ij=1,lmn2_size
       dij_up=half*(dijfock_vv_tmp(klmn_ij,1)+dijfock_vv_tmp(klmn_ij,4))
       dij_dn=half*(dijfock_vv_tmp(klmn_ij,1)-dijfock_vv_tmp(klmn_ij,4))
       dij_updn_r= half*dijfock_vv_tmp(klmn_ij,2)
       dij_updn_i=-half*dijfock_vv_tmp(klmn_ij,3)
       dijfock_vv(klmn1  ,1)= dij_up
       dijfock_vv(klmn1  ,2)= dij_dn
       dijfock_vv(klmn1  ,3)= dij_updn_r
       dijfock_vv(klmn1+1,3)= dij_updn_i
       dijfock_vv(klmn1  ,4)= dij_updn_r
       dijfock_vv(klmn1+1,4)=-dij_updn_i
       klmn1=klmn1+cplex_dij
     end do
   end if

 end do ! qphase

! ===== Core-valence contribution =====

 do ispden=1,pawrhoij%nsppol
   do iq=1,qphase
     iq0_dij=merge(0,cplex_dij*lmn2_size,iq==1)
     klmn1=iq0_dij+1
     do klmn_ij=1,lmn2_size
       dijfock_cv(klmn1,ispden)=pawtab%ex_cvij(klmn_ij)
       klmn1=klmn1+cplex_dij
     end do
   end do
 end do

!Antiferro case: copy up component into down one
 if (ndij==2.and.pawrhoij%nsppol==1) then
   dijfock_cv(:,2)=dijfock_cv(:,1)
 end if

!Apply mixing factors
 if (abs(hyb_mixing)>tol8) then
   dijfock_vv(:,:) = hyb_mixing*dijfock_vv(:,:)
 else if (abs(hyb_mixing_sr)>tol8) then
   dijfock_vv(:,:) = hyb_mixing_sr*dijfock_vv(:,:)
 end if
 dijfock_cv(:,:) = (hyb_mixing+hyb_mixing_sr)*dijfock_cv(:,:)

!Free temporary memory spaces
 LIBPAW_DEALLOCATE(dijfock_vv_tmp)

end subroutine pawdijfock
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawdijxc
!! NAME
!! pawdijxc
!!
!! FUNCTION
!! Compute the eXchange-Correlation contribution to the PAW pseudopotential strength Dij,
!! using densities and potential expressed on a (r,theta,phi) grid
!! (for one atom only):
!!   D_ij^XC= < Phi_i|Vxc( n1+ nc[+nhat])| Phi_j>
!!           -<tPhi_i|Vxc(tn1+tnc[+nhat])|tPhi_j>
!!           -Intg_omega [ Vxc(tn1+tnc[+nhat])(r). Sum_L(Qij^L(r)). dr]
!!
!! INPUTS
!!  cplex_dij=2 if dij is COMPLEX (as in the spin-orbit case), 1 if dij is REAL
!!  qphase=2 if dij contains a exp(-i.q.r) phase (as in the q<>0 RF case), 1 if not
!!  ndij= number of spin components
!!  nspden=number of spin density components
!!  nsppol=number of independent spin WF components
!!  pawang <type(pawang_type)>=paw angular mesh and related data, for current atom
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data, for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  usekden=1 if kinetic energy density contribution has to be included (mGGA)
!!  usexcnhat= 1 if compensation density is included in Vxc, 0 otherwise
!!  vxc1(qphase*mesh_size,angl_size,nspden)=all-electron on-site XC potential for current atom
!!                                   given on a (r,theta,phi) grid
!!  vxct1(qphase*mesh_size,angl_size,nspden)=all-electron on-site XC potential for current atom
!!                                    given on a (r,theta,phi) grid
!!  [vxctau1(qphase*mesh_size,angl_size,nspden)]=1st deriv. of XC energy wrt to kinetic energy density
!!                                               (all electron) - metaGGA only
!!  [vxcttau1(qphase*mesh_size,angl_size,nspden)]=1st deriv. of XC energy wrt to kinetic energy density
!!                                                (pseudo) - metaGGA only
!!
!! OUTPUT
!!  dijxc(cplex_dij*qphase*lmn2_size,ndij)=  D_ij^XC terms
!!    When Dij is complex (cplex_dij=2):
!!      dij(2*i-1,:) contains the real part, dij(2*i,:) contains the imaginary part
!!    When a exp(-i.q.r) phase is included (qphase=2):
!!      dij(1:cplex_dij*lmn2_size,:)
!!          contains the real part of the phase, i.e. D_ij*cos(q.r)
!!      dij(cplex_dij*lmn2_size+1:2*cplex_dij*lmn2_size,:)
!!          contains the imaginary part of the phase, i.e. D_ij*sin(q.r)
!!
!! PARENTS
!!      m_pawdij
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijxc(dijxc,cplex_dij,qphase,ndij,nspden,nsppol,&
&                   pawang,pawrad,pawtab,vxc1,vxct1,usexcnhat,usekden,&
&                   vxctau1,vxcttau1) ! optional

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex_dij,ndij,nspden,nsppol,qphase,usekden,usexcnhat
 type(pawang_type),intent(in) :: pawang
!arrays
 real(dp),intent(in) :: vxc1(:,:,:),vxct1(:,:,:)
 real(dp),intent(in),optional :: vxctau1(:,:,:),vxcttau1(:,:,:)
 real(dp),intent(out) :: dijxc(:,:)
 type(pawrad_type),intent(in) :: pawrad
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: angl_size,basis_size,idij,idijend
 integer :: ii,ij_size,ilm,ils,ils1,ilslm,ipts,ir,ir1,isel,ispden
 integer :: jlm,j0lm,klmn1,klmn2,klmn,klm,kln,l_size,lm0,lmax,lmin,lm_size,lmn2_size
 integer :: iln,jln,j0ln
 integer :: mesh_size,mm,nsploop
 real(dp) :: tmp,vi,vr,vxcijhat,vxcijhat_i,vxctauij
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: dijxc_idij(:),ff(:),gg(:)
 real(dp),allocatable :: vxcij1(:),vxcij2(:),vxctauij1(:),yylmr(:,:),yylmgr(:,:)

! *************************************************************************

!Useful data
 lm_size=pawtab%lcut_size**2
 lmn2_size=pawtab%lmn2_size
 basis_size=pawtab%basis_size
 ij_size=pawtab%ij_size
 l_size=pawtab%l_size
 mesh_size=pawtab%mesh_size
 angl_size=pawang%angl_size

!Check data consistency
 if (size(dijxc,1)/=cplex_dij*qphase*lmn2_size.or.size(dijxc,2)/=ndij) then
   msg='invalid sizes for Dijxc !'
   MSG_BUG(msg)
 end if
 if (size(vxc1,1)/=qphase*mesh_size.or.size(vxct1,1)/=qphase*mesh_size.or.&
&    size(vxc1,2)/=angl_size.or.size(vxct1,2)/=angl_size.or.&
&    size(vxc1,3)/=nspden.or.size(vxct1,3)/=nspden) then
   msg='invalid sizes for vxc1 or vxct1 !'
   MSG_BUG(msg)
 end if

!Check if MetaGGA is activated
 if (usekden==1) then
   if (.not.present(vxctau1)) then
     msg="vxctau1 needs to be present!"
     MSG_BUG(msg)
   else if (size(vxctau1)==0) then
     msg="vxctau1 needs to be allocated!"
     MSG_BUG(msg)
   end if
   if (.not.present(vxcttau1)) then
     msg="vxcttau1 needs to be present!"
     MSG_BUG(msg)
   else if (size(vxcttau1)==0) then
     msg="vxcttau1 needs to be allocated!"
     MSG_BUG(msg)
   end if
 end if

!Precompute products Ylm*Ylpmp (and Grad(Ylm).Grad(Ylpmp) if mGGA)
 lmax=1+maxval(pawtab%indklmn(4,1:lmn2_size))
 LIBPAW_ALLOCATE(yylmr,(lmax**2*(lmax**2+1)/2,angl_size))
 LIBPAW_ALLOCATE(yylmgr,(lmax**2*(lmax**2+1)/2,angl_size*usekden))
 do ipts=1,angl_size
   do jlm=1,lmax**2
     j0lm=jlm*(jlm-1)/2
     do ilm=1,jlm
       klm=j0lm+ilm
       yylmr(klm,ipts)=pawang%ylmr(ilm,ipts)*pawang%ylmr(jlm,ipts)
     end do
   end do
 end do
 if (usekden==1) then
   yylmgr(:,:)=zero
   do ipts=1,angl_size
     do jlm=1,lmax**2
       j0lm=jlm*(jlm-1)/2
       do ilm=1,jlm
         klm=j0lm+ilm
         do ii=1,3
           yylmgr(klm,ipts)=yylmgr(klm,ipts) &
&            +pawang%ylmrgr(ii,ilm,ipts)*pawang%ylmrgr(ii,jlm,ipts)
         end do
       end do
     end do
   end do
 end if

!Init memory
 dijxc=zero
 LIBPAW_ALLOCATE(dijxc_idij,(qphase*lmn2_size))
 LIBPAW_ALLOCATE(vxcij1,(qphase*ij_size))
 LIBPAW_ALLOCATE(vxcij2,(qphase*l_size))
 LIBPAW_ALLOCATE(vxctauij1,(qphase*ij_size*usekden))
 LIBPAW_ALLOCATE(ff,(mesh_size))
 LIBPAW_ALLOCATE(gg,(mesh_size))

!----------------------------------------------------------
!Loop over spin components
!----------------------------------------------------------
 nsploop=nsppol;if (ndij==4) nsploop=4
 do idij=1,nsploop
   if (idij<=nsppol.or.(nspden==4.and.idij<=3)) then

     idijend=idij+idij/3
     do ispden=idij,idijend

       dijxc_idij=zero

!      ----------------------------------------------------------
!      Loop on angular mesh
!      ----------------------------------------------------------
       do ipts=1,angl_size

!        ===== Vxc_ij_1 (tmp) =====
         vxcij1=zero
         if (qphase==1) then
           do kln=1,ij_size
             ff(1:mesh_size)= &
&               vxc1(1:mesh_size,ipts,ispden)*pawtab%phiphj(1:mesh_size,kln) &
&              -vxct1(1:mesh_size,ipts,ispden)*pawtab%tphitphj(1:mesh_size,kln)
             call simp_gen(vxcij1(kln),ff,pawrad)
           end do
           !if Meta GGA add 1/2*[<nabla_phi_i|vxctau1|nabla_phi_j>
           !                    -<nabla_tphi_i|vxcttau|nabla_tphi_j>]
           if (usekden==1) then
             do jln=1,basis_size
               j0ln=jln*(jln-1)/2
               do iln=1,jln
                 kln=j0ln+iln
                 ff(1:mesh_size)=vxctau1(1:mesh_size,ipts,ispden)*pawtab%phiphj(1:mesh_size,kln) &
&                               -vxcttau1(1:mesh_size,ipts,ispden)*pawtab%tphitphj(1:mesh_size,kln)
                 call simp_gen(vxctauij1(kln),ff,pawrad)
                 ff(1:mesh_size)=vxctau1(1:mesh_size,ipts,ispden) &
&                               *pawtab%nablaphi(1:mesh_size,iln)*pawtab%nablaphi(1:mesh_size,jln) &
&                               -vxcttau1(1:mesh_size,ipts,ispden) &
&                               *pawtab%tnablaphi(1:mesh_size,iln)*pawtab%tnablaphi(1:mesh_size,jln)
                 call simp_gen(vxctauij,ff,pawrad)
                 vxcij1(kln)=vxcij1(kln)+half*vxctauij
               end do
             end do
           end if

         else
           do kln=1,ij_size
             do ir=1,mesh_size
               ir1=2*ir
               ff(ir)= &
&                 vxc1(ir1-1,ipts,ispden)*pawtab%phiphj(ir,kln) &
&                -vxct1(ir1-1,ipts,ispden)*pawtab%tphitphj(ir,kln)
               gg(ir)= &
&                 vxc1(ir1,ipts,ispden)*pawtab%phiphj(ir,kln) &
&                -vxct1(ir1,ipts,ispden)*pawtab%tphitphj(ir,kln)
             end do
             call simp_gen(vxcij1(2*kln-1),ff,pawrad)
             call simp_gen(vxcij1(2*kln  ),gg,pawrad)
           end do
           !if Meta GGA add 1/2*[<nabla_phi_i|vxctau1|nabla_phi_j>
           !                    -<nabla_tphi_i|vxcttau|nabla_tphi_j>]
           if (usekden==1) then
             do jln=1,basis_size
               j0ln=jln*(jln-1)/2
               do iln=1,jln
                 kln=j0ln+iln
                 do ir=1,mesh_size
                   ir1=2*ir
                   ff(ir)=vxctau1(ir1-1,ipts,ispden)*pawtab%phiphj(ir,kln) &
&                        -vxcttau1(ir1-1,ipts,ispden)*pawtab%tphitphj(ir,kln)
                   gg(ir)=vxctau1(ir1,ipts,ispden)*pawtab%phiphj(ir,kln) &
&                        -vxcttau1(ir1,ipts,ispden)*pawtab%tphitphj(ir,kln)
                 end do
                 call simp_gen(vxctauij1(2*kln-1),ff,pawrad)
                 call simp_gen(vxctauij1(2*kln  ),gg,pawrad)
                 do ir=1,mesh_size
                   ir1=2*ir
                   ff(ir)=vxctau1(ir1-1,ipts,ispden) &
&                        *pawtab%nablaphi(ir,iln)*pawtab%nablaphi(ir,jln) &
&                        -vxcttau1(ir1-1,ipts,ispden) &
&                        *pawtab%tnablaphi(ir,iln)*pawtab%tnablaphi(ir,jln)
                   gg(ir)=vxctau1(ir1,ipts,ispden) &
&                        *pawtab%nablaphi(ir,iln)*pawtab%nablaphi(ir,jln) &
&                        -vxcttau1(ir1,ipts,ispden) &
&                        *pawtab%tnablaphi(ir,iln)*pawtab%tnablaphi(ir,jln)
                 end do
                 call simp_gen(vxctauij,ff,pawrad)
                 vxcij1(2*kln-1)=vxcij1(2*kln-1)+half*vxctauij
                 call simp_gen(vxctauij,gg,pawrad)
                 vxcij1(2*kln)=vxcij1(2*kln)+half*vxctauij
               end do
             end do
           end if

         end if

!        ===== Vxc_ij_2 (tmp) =====
         vxcij2=zero
         if (usexcnhat/=0) then
           if (qphase==1) then
             do ils=1,l_size
               ff(1:mesh_size)=vxct1(1:mesh_size,ipts,ispden) &
&                 *pawtab%shapefunc(1:mesh_size,ils) &
&                 *pawrad%rad(1:mesh_size)**2
               call simp_gen(vxcij2(ils),ff,pawrad)
             end do
           else
             do ils=1,l_size
               do ir=1,mesh_size
                 ir1=2*ir
                 tmp=pawtab%shapefunc(ir,ils)*pawrad%rad(ir)**2
                 ff(ir)=vxct1(ir1-1,ipts,ispden)*tmp
                 gg(ir)=vxct1(ir1  ,ipts,ispden)*tmp
               end do
               call simp_gen(vxcij2(2*ils-1),ff,pawrad)
               call simp_gen(vxcij2(2*ils  ),gg,pawrad)
             end do
           end if
         end if

!        ===== Integrate Vxc_ij_1 and Vxc_ij_2 over the angular mesh =====
!        ===== and accummulate in total Vxc_ij                       =====
         if (qphase==1) then
           do klmn=1,lmn2_size
             klm=pawtab%indklmn(1,klmn);kln=pawtab%indklmn(2,klmn)
             lmin=pawtab%indklmn(3,klmn);lmax=pawtab%indklmn(4,klmn)
             dijxc_idij(klmn)=dijxc_idij(klmn)+vxcij1(kln) &
&                            *pawang%angwgth(ipts)*yylmr(klm,ipts)*four_pi
             if (usekden==1) then
               dijxc_idij(klmn)=dijxc_idij(klmn)+half*vxctauij1(kln) &
&                              *pawang%angwgth(ipts)*yylmgr(klm,ipts)*four_pi
             end if
             if (usexcnhat/=0) then
               vxcijhat=zero
               do ils=lmin,lmax,2
                 lm0=ils**2+ils+1
                 vr=four_pi*pawang%angwgth(ipts)*vxcij2(ils+1)
                 do mm=-ils,ils
                   ilslm=lm0+mm;isel=pawang%gntselect(ilslm,klm)
                   if (isel>0) then
                     tmp=pawang%ylmr(ilslm,ipts)*pawtab%qijl(ilslm,klmn)
                     vxcijhat=vxcijhat+vr*tmp
                   end if
                 end do
               end do
               dijxc_idij(klmn)=dijxc_idij(klmn)-vxcijhat
             end if
           end do ! Loop klmn
         else
           klmn1=1
           do klmn=1,lmn2_size
             klm=pawtab%indklmn(1,klmn);kln=pawtab%indklmn(2,klmn)
             lmin=pawtab%indklmn(3,klmn);lmax=pawtab%indklmn(4,klmn)
             tmp=pawang%angwgth(ipts)*yylmr(klm,ipts)*four_pi
             dijxc_idij(klmn1  )=dijxc_idij(klmn1  )+vxcij1(2*kln-1)*tmp
             dijxc_idij(klmn1+1)=dijxc_idij(klmn1+1)+vxcij1(2*kln  )*tmp
             if (usekden==1) then
               tmp=pawang%angwgth(ipts)*yylmgr(klm,ipts)*four_pi
               dijxc_idij(klmn1  )=dijxc_idij(klmn1  )+half*vxctauij1(2*kln-1)*tmp
               dijxc_idij(klmn1+1)=dijxc_idij(klmn1+1)+half*vxctauij1(2*kln  )*tmp
             end if
             if (usexcnhat/=0) then
               vxcijhat=zero;vxcijhat_i=zero
               do ils=lmin,lmax,2
                 lm0=ils**2+ils+1;ils1=2*(ils+1)
                 vr=four_pi*pawang%angwgth(ipts)*vxcij2(ils1-1)
                 vi=four_pi*pawang%angwgth(ipts)*vxcij2(ils1  )
                 do mm=-ils,ils
                   ilslm=lm0+mm;isel=pawang%gntselect(ilslm,klm)
                   if (isel>0) then
                     tmp=pawang%ylmr(ilslm,ipts)*pawtab%qijl(ilslm,klmn)
                     vxcijhat  =vxcijhat  +vr*tmp
                     vxcijhat_i=vxcijhat_i+vi*tmp
                   end if
                 end do
               end do
               dijxc_idij(klmn1  )=dijxc_idij(klmn1  )-vxcijhat
               dijxc_idij(klmn1+1)=dijxc_idij(klmn1+1)-vxcijhat_i
             end if
             klmn1=klmn1+qphase
           end do ! Loop klmn
         end if

!      ----------------------------------------------------------
!      End loop on angular points
       end do

!      ----------------------------------------------------------
!      Deduce some part of Dij according to symmetries
!      ----------------------------------------------------------

       !if ispden=1 => real part of D^11_ij
       !if ispden=2 => real part of D^22_ij
       !if ispden=3 => real part of D^12_ij
       !if ispden=4 => imaginary part of D^12_ij
       klmn1=max(1,ispden-2);klmn2=1
       do klmn=1,lmn2_size
         dijxc(klmn1,idij)=dijxc_idij(klmn2)
         klmn1=klmn1+cplex_dij
         klmn2=klmn2+qphase
       end do
       if (qphase==2) then
         !Same storage with exp^(-i.q.r) phase
         klmn1=max(1,ispden-2)+lmn2_size*cplex_dij;klmn2=2
         do klmn=1,lmn2_size
           dijxc(klmn1,idij)=dijxc_idij(klmn2)
           klmn1=klmn1+cplex_dij
           klmn2=klmn2+qphase
         end do
       endif

     end do !ispden

   !Non-collinear: D_ij(:,4)=D^21_ij=D^12_ij^*
   else if (nspden==4.and.idij==4) then
     dijxc(:,idij)=dijxc(:,idij-1)
     if (cplex_dij==2) then
       do klmn=2,lmn2_size*cplex_dij,cplex_dij
         dijxc(klmn,idij)=-dijxc(klmn,idij)
       end do
       if (qphase==2) then
         do klmn=2+lmn2_size*cplex_dij,2*lmn2_size*cplex_dij,cplex_dij
           dijxc(klmn,idij)=-dijxc(klmn,idij)
         end do
       end if
     end if

   !Antiferro: D_ij(:,2)=D^down_ij=D^up_ij
   else if (nsppol==1.and.idij==2) then
     dijxc(:,idij)=dijxc(:,idij-1)
   end if

!----------------------------------------------------------
!End loop on spin density components
 end do

!Free temporary memory spaces
 LIBPAW_DEALLOCATE(yylmr)
 LIBPAW_DEALLOCATE(yylmgr)
 LIBPAW_DEALLOCATE(dijxc_idij)
 LIBPAW_DEALLOCATE(vxcij1)
 LIBPAW_DEALLOCATE(vxcij2)
 LIBPAW_DEALLOCATE(vxctauij1)
 LIBPAW_DEALLOCATE(ff)
 LIBPAW_DEALLOCATE(gg)

end subroutine pawdijxc
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawdijxcm
!! NAME
!! pawdijxcm
!!
!! FUNCTION
!! Compute the eXchange-Correlation contribution to the PAW pseudopotential strength Dij,
!! using densities and potential expressed as (l,m) spherical moments
!! (for one atom only):
!!   D_ij^XC= < Phi_i|Vxc( n1+ nc[+nhat])| Phi_j>
!!           -<tPhi_i|Vxc(tn1+tnc[+nhat])|tPhi_j>
!!           -Intg_omega [ Vxc(tn1+tnc[+nhat])(r). Sum_L(Qij^L(r)). dr]
!!
!! INPUTS
!!  cplex_dij=2 if dij is COMPLEX (as in the spin-orbit case), 1 if dij is REAL
!!  qphase=2 if dij contains a exp(-i.q.r) phase (as in the q<>0 RF case), 1 if not
!!  lmselect(lm_size)=select the non-zero LM-moments of on-site potentials
!!  ndij= number of spin components
!!  nspden=number of spin density components
!!  nsppol=number of independent spin WF components
!!  pawang <type(pawang_type)>=paw angular mesh and related data, for current atom
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data, for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  vxc1(qphase*mesh_size,lm_size,nspden)=all-electron on-site XC potential for current atom
!!                                 given on (l,m) spherical moments
!!  vxct1(qphase*mesh_size,lm_size,nspden)=all-electron on-site XC potential for current atom
!!                                  given on (l,m) spherical moments
!!  usexcnhat= 1 if compensation density is included in Vxc, 0 otherwise
!!
!! OUTPUT
!!  dijxc(cplex_dij*qphase*lmn2_size,ndij)=  D_ij^XC terms
!!    When Dij is complex (cplex_dij=2):
!!      dij(2*i-1,:) contains the real part, dij(2*i,:) contains the imaginary part
!!    When a exp(-i.q.r) phase is included (qphase=2):
!!      dij(1:cplex_dij*lmn2_size,:)
!!          contains the real part of the phase, i.e. D_ij*cos(q.r)
!!      dij(cplex_dij*lmn2_size+1:2*cplex_dij*lmn2_size,:)
!!          contains the imaginary part of the phase, i.e. D_ij*sin(q.r)
!!
!! PARENTS
!!      m_pawdij
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijxcm(dijxc,cplex_dij,qphase,lmselect,ndij,nspden,nsppol,&
&                    pawang,pawrad,pawtab,vxc1,vxct1,usexcnhat)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex_dij,ndij,nspden,nsppol,qphase,usexcnhat
 type(pawang_type),intent(in) :: pawang
!arrays
 logical :: lmselect(:)
 real(dp),intent(in) :: vxc1(:,:,:),vxct1(:,:,:)
 real(dp),intent(out) :: dijxc(:,:)
 type(pawrad_type),intent(in) :: pawrad
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: idij,idijend,ij_size,ir,ir1,isel,ispden,klm,klm1,klmn,klmn1,klmn2,kln
 integer :: lm_size,lmn2_size,ll,mesh_size,nsploop
 real(dp) :: tmp,vxcij2,vxcij2_i
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: dijxc_idij(:),ff(:),gg(:),vxcij1(:)

! *************************************************************************

!Useful data
 lm_size=pawtab%lcut_size**2
 lmn2_size=pawtab%lmn2_size
 ij_size=pawtab%ij_size
 mesh_size=pawtab%mesh_size

!Check data consistency
 if (size(dijxc,1)/=cplex_dij*qphase*lmn2_size.or.size(dijxc,2)/=ndij) then
   msg='invalid sizes for Dijxc !'
   MSG_BUG(msg)
 end if
 if (size(lmselect)/=lm_size) then
   msg='invalid size for lmselect !'
   MSG_BUG(msg)
 end if
 if (size(vxc1,1)/=qphase*mesh_size.or.size(vxct1,1)/=qphase*mesh_size.or.&
&    size(vxc1,2)/=lm_size.or.size(vxct1,2)/=lm_size.or.&
&    size(vxc1,3)/=nspden.or.size(vxct1,3)/=nspden) then
   msg='invalid sizes for vxc1 or vxct1 !'
   MSG_BUG(msg)
 end if

!Init memory
 dijxc=zero
 LIBPAW_ALLOCATE(dijxc_idij,(qphase*lmn2_size))
 LIBPAW_ALLOCATE(vxcij1,(qphase*ij_size))
 LIBPAW_ALLOCATE(ff,(mesh_size))
 LIBPAW_ALLOCATE(gg,(mesh_size))

!----------------------------------------------------------
!Loop over spin components
!----------------------------------------------------------
 nsploop=nsppol;if (ndij==4) nsploop=4
 do idij=1,nsploop
   if (idij<=nsppol.or.(nspden==4.and.idij<=3)) then

     idijend=idij+idij/3
     do ispden=idij,idijend

       dijxc_idij=zero

!      ----------------------------------------------------------
!      Summing over (l,m) moments
!      ----------------------------------------------------------
       do klm=1,lm_size
         if (lmselect(klm)) then

!          ===== Vxc_ij_1 (tmp) =====
           vxcij1=zero
           if (qphase==1) then
             do kln=1,ij_size
               ff(1:mesh_size)= &
&                 vxc1(1:mesh_size,klm,ispden)*pawtab%phiphj(1:mesh_size,kln) &
&                -vxct1(1:mesh_size,klm,ispden)*pawtab%tphitphj(1:mesh_size,kln)
               call simp_gen(vxcij1(kln),ff,pawrad)
             end do
           else ! qphase==2
             do kln=1,ij_size
               do ir=1,mesh_size
                 ir1=2*ir
                 ff(ir)= &
&                   vxc1(ir1-1,klm,ispden)*pawtab%phiphj(ir,kln) &
&                  -vxct1(ir1-1,klm,ispden)*pawtab%tphitphj(ir,kln)
                 gg(ir)= &
&                   vxc1(ir1,klm,ispden)*pawtab%phiphj(ir,kln) &
&                  -vxct1(ir1,klm,ispden)*pawtab%tphitphj(ir,kln)
               end do
               call simp_gen(vxcij1(2*kln-1),ff,pawrad)
               call simp_gen(vxcij1(2*kln  ),gg,pawrad)
             end do
           end if

!          ===== Vxc_ij_2 (tmp) =====
           vxcij2=zero;vxcij2_i=zero
           if (usexcnhat/=0) then
             ll=1+int(sqrt(dble(klm)-0.1_dp))
             if (qphase==1) then
               ff(1:mesh_size)=vxct1(1:mesh_size,klm,ispden) &
&                             *pawtab%shapefunc(1:mesh_size,ll) &
&                             *pawrad%rad(1:mesh_size)**2
               call simp_gen(vxcij2,ff,pawrad)
             else ! qphase==2
               do ir=1,mesh_size
                 ir1=2*ir
                 tmp=pawtab%shapefunc(ir,ll)*pawrad%rad(ir)**2
                 ff(ir)=vxct1(ir1-1,klm,ispden)*tmp
                 gg(ir)=vxct1(ir1  ,klm,ispden)*tmp
               end do
               call simp_gen(vxcij2  ,ff,pawrad)
               call simp_gen(vxcij2_i,gg,pawrad)
             end if
           end if

!          ===== Accumulate over klm moments Vxc_ij_1 and Vxc_ij_2 =====
!          ===== into total Vxc_ij                                 =====
           if (qphase==1) then
             do klmn=1,lmn2_size
               klm1=pawtab%indklmn(1,klmn)
               kln=pawtab%indklmn(2,klmn)
               isel=pawang%gntselect(klm,klm1)
               if (isel>0) &
&                dijxc_idij(klmn)=dijxc_idij(klmn)+vxcij1(kln)*pawang%realgnt(isel)
               if (usexcnhat/=0) &
                 dijxc_idij(klmn)=dijxc_idij(klmn)-pawtab%qijl(klm,klmn)*vxcij2
             end do ! Loop klmn
           else ! qphase==2
             klmn1=1
             do klmn=1,lmn2_size
               klm1=pawtab%indklmn(1,klmn)
               kln=pawtab%indklmn(2,klmn)
               isel=pawang%gntselect(klm,klm1)
               if (isel>0) then
                 dijxc_idij(klmn1  )=dijxc_idij(klmn1) &
&                                   +vxcij1(2*kln-1)*pawang%realgnt(isel)
                 dijxc_idij(klmn1+1)=dijxc_idij(klmn1+1) &
&                                   +vxcij1(2*kln  )*pawang%realgnt(isel)
               end if
               if (usexcnhat/=0) then
                 dijxc_idij(klmn1  )=dijxc_idij(klmn1) &
&                                   -pawtab%qijl(klm,klmn)*vxcij2
                 dijxc_idij(klmn1+1)=dijxc_idij(klmn1+1) &
&                                   -pawtab%qijl(klm,klmn)*vxcij2_i
               end if
               klmn1=klmn1+qphase
             end do ! Loop klmn
           end if

         end if ! klm selection
       end do  ! Loop klm

!      ----------------------------------------------------------
!      Deduce some part of Dij according to symmetries
!      ----------------------------------------------------------

       !if ispden=1 => real part of D^11_ij
       !if ispden=2 => real part of D^22_ij
       !if ispden=3 => real part of D^12_ij
       !if ispden=4 => imaginary part of D^12_ij
       klmn1=max(1,ispden-2);klmn2=1
       do klmn=1,lmn2_size
         dijxc(klmn1,idij)=dijxc_idij(klmn2)
         klmn1=klmn1+cplex_dij
         klmn2=klmn2+qphase
       end do
       if (qphase==2) then
         !Same storage with exp^(-i.q.r) phase
         klmn1=max(1,ispden-2)+lmn2_size*cplex_dij;klmn2=2
         do klmn=1,lmn2_size
           dijxc(klmn1,idij)=dijxc_idij(klmn2)
           klmn1=klmn1+cplex_dij
           klmn2=klmn2+qphase
         end do
       endif

     end do !ispden

   !Non-collinear: D_ij(:,4)=D^21_ij=D^12_ij^*
   else if (nspden==4.and.idij==4) then
     dijxc(:,idij)=dijxc(:,idij-1)
     if (cplex_dij==2) then
       do klmn=2,lmn2_size*cplex_dij,cplex_dij
         dijxc(klmn,idij)=-dijxc(klmn,idij)
       end do
       if (qphase==2) then
         do klmn=2+lmn2_size*cplex_dij,2*lmn2_size*cplex_dij,cplex_dij
           dijxc(klmn,idij)=-dijxc(klmn,idij)
         end do
       end if
     end if

   !Antiferro: D_ij(:,2)=D^down_ij=D^up_ij
   else if (nsppol==1.and.idij==2) then
     dijxc(:,idij)=dijxc(:,idij-1)
   end if

!----------------------------------------------------------
!End loop on spin density components
 end do

!Free temporary memory spaces
 LIBPAW_DEALLOCATE(dijxc_idij)
 LIBPAW_DEALLOCATE(vxcij1)
 LIBPAW_DEALLOCATE(ff)
 LIBPAW_DEALLOCATE(gg)

end subroutine pawdijxcm
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawdijhat
!! NAME
!! pawdijhat
!!
!! FUNCTION
!! Compute the "hat" contribution to the PAW pseudopotential strength Dij,
!! i.e. the compensation charge contribution (for one atom only):
!!   D_ij^hat=Intg_R [ V(r). Sum_L(Qij^L(r)). dr]
!!
!! INPUTS
!!  cplex_dij=2 if dij is COMPLEX (as in the spin-orbit case), 1 if dij is REAL
!!  qphase=2 if dij contains a exp(-i.q.r) phase (as in the q<>0 RF case), 1 if not
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!  iatom=absolute index of current atom (between 1 and natom)
!!  natom=total number of atoms
!!  ndij= number of spin components
!!  ngrid=number of points of the real space grid (FFT, WVL, ...) treated by current proc
!!  ngridtot=total number of points of the real space grid (FFT, WVL, ...)
!!           For the FFT grid, thi should be equal to ngfft1*ngfft2*ngfft3
!!  nspden=number of spin density components
!!  nsppol=number of independent spin WF components
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab<type(pawfgrtab_type)>=atomic data given on fine rectangular grid for current atom
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data, for current atom
!!  Pot(qphase*ngrid,nspden)=potential on real space grid
!!  qphon(3)=(RF calculations only) - wavevector of the phonon
!!  ucvol=unit cell volume
!!  xred(3,my_natom)= reduced atomic coordinates
!!
!! OUTPUT
!!  dijhat(cplex_dij*qphase*lmn2_size,ndij)= D_ij^hat terms
!!    When Dij is complex (cplex_dij=2):
!!      dij(2*i-1,:) contains the real part, dij(2*i,:) contains the imaginary part
!!    When a exp(-i.q.r) phase is included (qphase=2):
!!      dij(1:cplex_dij*lmn2_size,:)
!!          contains the real part of the phase, i.e. D_ij*cos(q.r)
!!      dij(cplex_dij*lmn2_size+1:2*cplex_dij*lmn2_size,:)
!!          contains the imaginary part of the phase, i.e. D_ij*sin(q.r)
!!
!! PARENTS
!!      fock_getghc,m_pawdij
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijhat(dijhat,cplex_dij,qphase,gprimd,iatom,&
&                    natom,ndij,ngrid,ngridtot,nspden,nsppol,pawang,pawfgrtab,&
&                    pawtab,Pot,qphon,ucvol,xred,&
&                    mpi_comm_grid) ! Optional argument

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex_dij,iatom,natom,ndij
 integer,intent(in) :: ngrid,ngridtot,nspden,nsppol,qphase
 integer,intent(in),optional :: mpi_comm_grid
 real(dp),intent(in) :: ucvol
 type(pawang_type),intent(in) :: pawang
 type(pawfgrtab_type),intent(inout) :: pawfgrtab
!arrays
 real(dp),intent(in) :: gprimd(3,3),Pot(qphase*ngrid,nspden),qphon(3),xred(3,natom)
 real(dp),intent(out) :: dijhat(:,:)
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: ic,idij,idijend,ier,ils,ilslm,ilslm1,isel,ispden,jc,klm,klmn,klmn1,klmn2
 integer :: lm0,lm_size,lmax,lmin,lmn2_size,mm,my_comm_grid,nfgd,nsploop,optgr0
 logical :: has_qphase,qne0
 real(dp) :: vi,vr
 character(len=500) :: msg
!arrays
 real(dp) :: rdum1(1),rdum2(2)
 real(dp),allocatable :: dijhat_idij(:),prod(:)

! *************************************************************************

!Useful data
 lm_size=pawtab%lcut_size**2
 lmn2_size=pawtab%lmn2_size
 nfgd=pawfgrtab%nfgd
 qne0=(qphon(1)**2+qphon(2)**2+qphon(3)**2>=1.d-15)
 has_qphase=(qne0.and.qphase==2)
 my_comm_grid=xmpi_comm_self;if (present(mpi_comm_grid)) my_comm_grid=mpi_comm_grid

!Check data consistency
 if (size(dijhat,1)/=cplex_dij*qphase*lmn2_size.or.size(dijhat,2)/=ndij) then
   msg='invalid sizes for Dijhat !'
   MSG_BUG(msg)
 end if

!Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
 if (pawfgrtab%gylm_allocated==0) then
   if (allocated(pawfgrtab%gylm))  then
     LIBPAW_DEALLOCATE(pawfgrtab%gylm)
   end if
   LIBPAW_ALLOCATE(pawfgrtab%gylm,(nfgd,lm_size))
   pawfgrtab%gylm_allocated=2;optgr0=1
   call pawgylm(pawfgrtab%gylm,rdum1,rdum2,lm_size,nfgd,optgr0,0,0,pawtab,pawfgrtab%rfgd)
 end if

!Eventually compute exp(i.q.r) factors for the current atom (if not already done)
 if (has_qphase.and.pawfgrtab%expiqr_allocated==0) then
   if (pawfgrtab%rfgd_allocated==0) then
     msg='pawfgrtab()%rfgd array must be allocated  !'
     MSG_BUG(msg)
   end if
   if (allocated(pawfgrtab%expiqr))  then
     LIBPAW_DEALLOCATE(pawfgrtab%expiqr)
   end if
   LIBPAW_ALLOCATE(pawfgrtab%expiqr,(2,nfgd))
   call pawexpiqr(pawfgrtab%expiqr,gprimd,nfgd,qphon,pawfgrtab%rfgd,xred(:,iatom))
   pawfgrtab%expiqr_allocated=2
 end if

!Init memory
 dijhat=zero
 LIBPAW_ALLOCATE(prod,(qphase*lm_size))
 LIBPAW_ALLOCATE(dijhat_idij,(qphase*lmn2_size))

!----------------------------------------------------------
!Loop over spin components
!----------------------------------------------------------
 nsploop=nsppol;if (ndij==4) nsploop=4
 do idij=1,nsploop
   if (idij<=nsppol.or.(nspden==4.and.idij<=3)) then

     idijend=idij+idij/3
     do ispden=idij,idijend

!      ------------------------------------------------------
!      Compute Int[V(r).g_l(r).Y_lm(r)]
!      ------------------------------------------------------
!       Note for non-collinear magnetism:
!          We compute Int[V^(alpha,beta)(r).g_l(r).Y_lm(r)]
!          Remember: if nspden=4, V is stored as : V^11, V^22, V^12, i.V^21

       prod=zero

!      ===== Standard case ============================
       if (.not.has_qphase) then
         if (qphase==1) then
           do ilslm=1,lm_size
             do ic=1,nfgd
               vr=Pot(pawfgrtab%ifftsph(ic),ispden)
               prod(ilslm)=prod(ilslm)+vr*pawfgrtab%gylm(ic,ilslm)
             end do
           end do
         else
           ilslm1=1
           do ilslm=1,lm_size
             do ic=1,nfgd
               jc=2*pawfgrtab%ifftsph(ic)
               vr=Pot(jc-1,ispden);vi=Pot(jc,ispden)
               prod(ilslm1  )=prod(ilslm1  )+vr*pawfgrtab%gylm(ic,ilslm)
               prod(ilslm1+1)=prod(ilslm1+1)+vi*pawfgrtab%gylm(ic,ilslm)
             end do
             ilslm1=ilslm1+qphase
           end do
         end if

!      ===== Including Exp(iqr) phase (DFPT only) =====
       else
         if (qphase==1) then
           do ilslm=1,lm_size
             do ic=1,nfgd
               vr=Pot(pawfgrtab%ifftsph(ic),ispden)
               prod(ilslm)=prod(ilslm)+vr*pawfgrtab%gylm(ic,ilslm)&
&                                        *pawfgrtab%expiqr(1,ic)
             end do
           end do
         else
           ilslm1=1
           do ilslm=1,lm_size
             do ic=1,nfgd
               jc=2*pawfgrtab%ifftsph(ic)
               vr=Pot(jc-1,ispden);vi=Pot(jc,ispden)
               prod(ilslm1  )=prod(ilslm1  )+pawfgrtab%gylm(ic,ilslm)&
&                *(vr*pawfgrtab%expiqr(1,ic)-vi*pawfgrtab%expiqr(2,ic))
               prod(ilslm1+1)=prod(ilslm1+1)+pawfgrtab%gylm(ic,ilslm)&
&                *(vr*pawfgrtab%expiqr(2,ic)+vi*pawfgrtab%expiqr(1,ic))
             end do
             ilslm1=ilslm1+qphase
           end do
         end if
       end if

!      Scaling factor (unit volume)
       prod=prod*ucvol/dble(ngridtot)

!      Reduction in case of parallelism
       if (xmpi_comm_size(my_comm_grid)>1) then
         call xmpi_sum(prod,my_comm_grid,ier)
       end if

!      ----------------------------------------------------------
!      Compute Sum_(i,j)_LM { q_ij^L Int[V(r).g_l(r).Y_lm(r)] }
!      ----------------------------------------------------------
!        Note for non-collinear magnetism:
!          We compute Sum_(i,j)_LM { q_ij^L Int[V^(alpha,beta)(r).g_l(r).Y_lm(r)] }

       dijhat_idij=zero

       if (qphase==1) then
         do klmn=1,lmn2_size
           klm =pawtab%indklmn(1,klmn)
           lmin=pawtab%indklmn(3,klmn)
           lmax=pawtab%indklmn(4,klmn)
           do ils=lmin,lmax,2
             lm0=ils**2+ils+1
             do mm=-ils,ils
               ilslm=lm0+mm;isel=pawang%gntselect(lm0+mm,klm)
               if (isel>0) dijhat_idij(klmn)=dijhat_idij(klmn) &
&                  +prod(ilslm)*pawtab%qijl(ilslm,klmn)
             end do
           end do
         end do
       else
         do klmn=1,lmn2_size
           klmn1=2*klmn-1
           klm =pawtab%indklmn(1,klmn)
           lmin=pawtab%indklmn(3,klmn)
           lmax=pawtab%indklmn(4,klmn)
           do ils=lmin,lmax,2
             lm0=ils**2+ils+1
             do mm=-ils,ils
               ilslm=lm0+mm;ilslm1=2*ilslm;isel=pawang%gntselect(lm0+mm,klm)
               if (isel>0) dijhat_idij(klmn1:klmn1+1)=dijhat_idij(klmn1:klmn1+1) &
&                  +prod(ilslm1-1:ilslm1)*pawtab%qijl(ilslm,klmn)
             end do
           end do
         end do
       end if

!      ----------------------------------------------------------
!      Deduce some part of Dij according to symmetries
!      ----------------------------------------------------------

       !if ispden=1 => real part of D^11_ij
       !if ispden=2 => real part of D^22_ij
       !if ispden=3 => real part of D^12_ij
       !if ispden=4 => imaginary part of D^12_ij
       klmn1=max(1,ispden-2);klmn2=1
       do klmn=1,lmn2_size
         dijhat(klmn1,idij)=dijhat_idij(klmn2)
         klmn1=klmn1+cplex_dij
         klmn2=klmn2+qphase
       end do
       if (qphase==2) then
         !Same storage with exp^(-i.q.r) phase
         klmn1=max(1,ispden-2)+lmn2_size*cplex_dij;klmn2=2
         do klmn=1,lmn2_size
           dijhat(klmn1,idij)=dijhat_idij(klmn2)
           klmn1=klmn1+cplex_dij
           klmn2=klmn2+qphase
         end do
       endif

     end do !ispden

   !Non-collinear: D_ij(:,4)=D^21_ij=D^12_ij^*
   else if (nspden==4.and.idij==4) then
     dijhat(:,idij)=dijhat(:,idij-1)
     if (cplex_dij==2) then
       do klmn=2,lmn2_size*cplex_dij,cplex_dij
         dijhat(klmn,idij)=-dijhat(klmn,idij)
       end do
       if (qphase==2) then
         do klmn=2+lmn2_size*cplex_dij,2*lmn2_size*cplex_dij,cplex_dij
           dijhat(klmn,idij)=-dijhat(klmn,idij)
         end do
       end if
     end if

   !Antiferro: D_ij(:,2)=D^down_ij=D^up_ij
   else if (nsppol==1.and.idij==2) then
     dijhat(:,idij)=dijhat(:,idij-1)
   end if

!----------------------------------------------------------
!End loop on spin density components
 end do

!Free temporary memory spaces
 LIBPAW_DEALLOCATE(prod)
 LIBPAW_DEALLOCATE(dijhat_idij)
 if (pawfgrtab%gylm_allocated==2) then
   LIBPAW_DEALLOCATE(pawfgrtab%gylm)
   LIBPAW_ALLOCATE(pawfgrtab%gylm,(0,0))
   pawfgrtab%gylm_allocated=0
 end if
 if (pawfgrtab%expiqr_allocated==2) then
   LIBPAW_DEALLOCATE(pawfgrtab%expiqr)
   LIBPAW_ALLOCATE(pawfgrtab%expiqr,(0,0))
   pawfgrtab%expiqr_allocated=0
 end if

end subroutine pawdijhat
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawdijnd
!! NAME
!! pawdijnd
!!
!! FUNCTION
!! Compute the nuclear dipole contribution to the PAW
!! pseudopotential strength Dij
!! (for one atom only)
!!
!! INPUTS
!!  cplex_dij=2 if dij is COMPLEX (as in the spin-orbit case), 1 if dij is REAL
!!  ndij= number of spin components
!!  nucdipmom(3) nuclear magnetic dipole moment for current atom
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data, for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data, for current atom
!!
!! OUTPUT
!!  dijnd(cplex_dij*lmn2_size,ndij)= nuclear dipole moment Dij terms
!!      cplex_dij=2 must be 2
!!      dij(2*i-1,:) contains the real part,
!!      dij(2*i  ,:) contains the imaginary part
!!
!! NOTES
!!   On-site contribution of a nuclear magnetic dipole moment at $R$. Hamiltonian is
!!   $H=(1/2m_e)(p - q_e A)^2 + V$ in SI units, and vector potential $A$ is
!!   $A=(\mu_0/4\pi) m\times (r-R)/|r-R|^3 = (\mu_0/4\pi) L_R\cdot m/|r-R|^3$ where
!!   $L_R$ is the on-site orbital angular momentum and $m$ is the nuclear magnetic
!!   dipole moment. Second order term in A is ignored. In atomic units the on-site term
!!   is \alpha^2 L_R\cdot m/|r-R|^3, where \alpha is the fine structure constant.
!!
!!
!! PARENTS
!!      m_pawdij,pawdenpot
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijnd(dijnd,cplex_dij,ndij,nucdipmom,pawrad,pawtab)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex_dij,ndij
 type(pawrad_type),intent(in) :: pawrad
 type(pawtab_type),target,intent(in) :: pawtab
!arrays
 real(dp),intent(out) :: dijnd(:,:)
 real(dp),intent(in) :: nucdipmom(3)

!Local variables ---------------------------------------
!scalars
 integer :: idir,ilmn,il,im,iln,ilm,jlmn,jl,jm,jlm,jln,klmn,kln,mesh_size
 real(dp) :: intgr3
 complex(dpc) :: lms
 logical :: ndmom
!arrays
 integer, LIBPAW_CONTIGUOUS pointer :: indlmn(:,:)
 real(dp),allocatable :: ff(:)
 character(len=500) :: msg

! *************************************************************************

!Useful data
 indlmn => pawtab%indlmn
 mesh_size=pawtab%mesh_size
 LIBPAW_ALLOCATE(ff,(mesh_size))

!Check data consistency
 if (cplex_dij/=2) then
   msg='cplex_dij must be 2 for nuclear dipole moments !'
   MSG_BUG(msg)
 end if
 if (size(dijnd,1)/=cplex_dij*pawtab%lmn2_size.or.size(dijnd,2)/=ndij) then
   msg='invalid sizes for Dijnd !'
   MSG_BUG(msg)
 end if

 dijnd = zero
 ndmom=(any(abs(nucdipmom)>tol8))

 if (ndmom) then ! only do the computation if at least one component of nuclear dipole is nonzero

!  loop over basis state pairs for this type
   do jlmn=1,pawtab%lmn_size
     jl=indlmn(1,jlmn)
     jm=indlmn(2,jlmn)
     jlm=indlmn(4,jlmn)
     jln=indlmn(5,jlmn)
     do ilmn=1,pawtab%lmn_size
       il=indlmn(1,ilmn)
       im=indlmn(2,ilmn)
       iln=indlmn(5,ilmn)
       ilm=indlmn(4,ilmn)
       klmn=max(jlmn,ilmn)*(max(jlmn,ilmn)-1)/2 + min(jlmn,ilmn)
       kln = pawtab%indklmn(2,klmn)

  !    Computation of (<phi_i|phi_j>-<tphi_i|tphi_j>)/r^3 radial integral

       ff(2:mesh_size)=(pawtab%phiphj(2:mesh_size,kln)-&
  &     pawtab%tphitphj(2:mesh_size,kln))/pawrad%rad(2:mesh_size)**3
       call pawrad_deducer0(ff,mesh_size,pawrad)
       call simp_gen(intgr3,ff,pawrad)

       do idir = 1, 3

! matrix element <S il im|L_idir|S jl jm>
         call slxyzs(il,im,idir,jl,jm,lms)

         dijnd(2*klmn-1,1) = dijnd(2*klmn-1,1) + &
              & intgr3*dreal(lms)*nucdipmom(idir)*FineStructureConstant2
         dijnd(2*klmn,1) = dijnd(2*klmn,1) + &
              & intgr3*dimag(lms)*nucdipmom(idir)*FineStructureConstant2

       end do

     end do ! end loop over ilmn
   end do ! end loop over jlmn

! in case of ndij > 1, note that there is no spin-flip in this term
! so therefore down-down = up-up, and up-down and down-up terms are still zero
   if(ndij > 1) dijnd(:,2)=dijnd(:,1)

 end if ! end check for a nonzero nuclear dipole moment

 LIBPAW_DEALLOCATE(ff)

end subroutine pawdijnd
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawdijso
!! NAME
!! pawdijso
!!
!! FUNCTION
!! Compute the spin-orbit contribution to the PAW
!! pseudopotential strength Dij
!! (for one atom only)
!!
!! INPUTS
!!  cplex_dij=2 if dij is COMPLEX (as in the spin-orbit case), 1 if dij is REAL
!!  qphase=2 if dij contains a exp(-i.q.r) phase (as in the q<>0 RF case), 1 if not
!!  ndij= number of spin components for Dij^SO
!!  nspden=number of spin density components
!!  paw_an <type(paw_an_type)>=paw arrays given on angular mesh, for current atom
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data, for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data, for current atom
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  spnorbscl=scaling factor for spin-orbit coupling
!!  vh1(qphase*mesh_size,v_size,nspden)=all-electron on-site Hartree potential for current atom
!!                     only spherical moment is used
!!  vxc1(qphase*mesh_size,v_size,nspden)=all-electron on-site XC potential for current atom
!!                                given on a (r,theta,phi) grid (v_size=angl_size)
!!                                or on (l,m) spherical moments (v_size=lm_size)
!!
!! OUTPUT
!!  dijso(cplex_dij*qphase*lmn2_size,ndij)= spin-orbit Dij terms
!!    Dij^SO is complex, so cplex_dij=2 must be 2:
!!      dij(2*i-1,:) contains the real part
!!      dij(2*i,:) contains the imaginary part
!!    Dij^SO is represented with 4 components:
!!      dijso(:,:,1) contains Dij_SO^up-up
!!      dijso(:,:,2) contains Dij_SO^dn-dn
!!      dijso(:,:,3) contains Dij_SO^up-dn
!!      dijso(:,:,4) contains Dij_SO^dn-up
!!    When a exp(-i.q.r) phase is included (qphase=2):
!!      dij(1:cplex_dij*lmn2_size,:)
!!          contains the real part of the phase, i.e. D_ij*cos(q.r)
!!      dij(cplex_dij*lmn2_size+1:2*cplex_dij*lmn2_size,:)
!!          contains the imaginary part of the phase, i.e. D_ij*sin(q.r)
!!
!! PARENTS
!!      m_pawdij,pawdenpot
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijso(dijso,cplex_dij,qphase,ndij,nspden,&
&                   pawang,pawrad,pawtab,pawxcdev,spnorbscl,vh1,vxc1)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex_dij,ndij,nspden,pawxcdev,qphase
 real(dp), intent(in) :: spnorbscl
 type(pawang_type),intent(in) :: pawang
!arrays
 real(dp),intent(out) :: dijso(:,:)
 real(dp),intent(in) :: vh1(:,:,:),vxc1(:,:,:)
 type(pawrad_type),intent(in) :: pawrad
 type(pawtab_type),target,intent(in) :: pawtab
!Local variables ---------------------------------------
!scalars
 integer :: angl_size,idij,ij_size,ilm,ipts,ispden,jlm,klm,klmn,klmn1,kln
 integer :: lm_size,lmn2_size,mesh_size,nsploop
 real(dp), parameter :: HalfFineStruct2=half/InvFineStruct**2
 real(dp) :: fact
 character(len=500) :: msg
!arrays
 integer, pointer :: indklmn(:,:)
 real(dp),allocatable :: dijso_rad(:),dv1dr(:),ff(:)

! *************************************************************************

!Useful data
 lm_size=pawtab%lcut_size**2
 lmn2_size=pawtab%lmn2_size
 ij_size=pawtab%ij_size
 angl_size=pawang%angl_size
 mesh_size=pawtab%mesh_size
 indklmn => pawtab%indklmn
 nsploop=4

!Check data consistency
 if (qphase/=1) then
   msg='qphase=2 not yet available in pawdijso!'
   MSG_BUG(msg)
 end if
 if (cplex_dij/=2) then
   msg='cplex_dij must be 2 for spin-orbit coupling!'
   MSG_BUG(msg)
 end if
 if (ndij/=4) then
   msg='ndij must be 4 for spin-orbit coupling!'
   MSG_BUG(msg)
 end if
 if (pawang%use_ls_ylm==0) then
   msg='pawang%use_ls_ylm should be /=0!'
   MSG_BUG(msg)
 end if
 if (size(dijso,1)/=cplex_dij*qphase*lmn2_size.or.size(dijso,2)/=ndij) then
   msg='invalid sizes for DijSO!'
   MSG_BUG(msg)
 end if
 if (size(vh1,1)/=qphase*mesh_size.or.size(vh1,2)<1.or.size(vh1,3)<1) then
   msg='invalid sizes for vh1!'
   MSG_BUG(msg)
 end if
 if (size(vxc1,1)/=qphase*mesh_size.or.size(vxc1,3)/=nspden.or.&
&   (size(vxc1,2)/=angl_size.and.pawxcdev==0).or.&
&   (size(vxc1,2)/=lm_size.and.pawxcdev/=0)) then
   msg='invalid sizes for vxc1!'
   MSG_BUG(msg)
 end if

!------------------------------------------------------------------------
!----------- Allocations and initializations
!------------------------------------------------------------------------

!Eventually compute <Phi_i|1/r.dV/dr|Phi_j>*alpha2/2*Y_00 (for spin-orbit)
 LIBPAW_ALLOCATE(dv1dr,(mesh_size))
 LIBPAW_ALLOCATE(dijso_rad,(ij_size))
 LIBPAW_ALLOCATE(ff,(mesh_size))
 fact=one/sqrt(four_pi) ! Y_00
 if (pawxcdev/=0) then
   if (nspden==1) then
     ff(1:mesh_size)=vxc1(1:mesh_size,1,1)
   else
     ff(1:mesh_size)=half*(vxc1(1:mesh_size,1,1)+vxc1(1:mesh_size,1,2))
   end if
 else
   ff(1:mesh_size)=zero
   if (nspden==1) then
     do ipts=1,angl_size
       ff(1:mesh_size)=ff(1:mesh_size) &
&          +vxc1(1:mesh_size,ipts,1)*pawang%angwgth(ipts)
     end do
   else
     do ipts=1,angl_size
       ff(1:mesh_size)=ff(1:mesh_size) &
&       +half*(vxc1(1:mesh_size,ipts,1)+vxc1(1:mesh_size,ipts,2)) &
&       *pawang%angwgth(ipts)
     end do
   end if
   ff(1:mesh_size)=sqrt(four_pi)*ff(1:mesh_size)
 end if
 ff(1:mesh_size)=fact*(ff(1:mesh_size)+vh1(1:mesh_size,1,1))
 call nderiv_gen(dv1dr,ff,pawrad)
 dv1dr(2:mesh_size)=HalfFineStruct2*(one/(one-ff(2:mesh_size)/InvFineStruct**2)) &
& *dv1dr(2:mesh_size)/pawrad%rad(2:mesh_size)
 call pawrad_deducer0(dv1dr,mesh_size,pawrad)
 do kln=1,ij_size
   ff(1:mesh_size)= dv1dr(1:mesh_size)*pawtab%phiphj(1:mesh_size,kln)
   call simp_gen(dijso_rad(kln),ff,pawrad)
 end do
 LIBPAW_DEALLOCATE(dv1dr)
 LIBPAW_DEALLOCATE(ff)
 dijso_rad(:)=spnorbscl*dijso_rad(:)

!------------------------------------------------------------------------
!----- Loop over density components
!------------------------------------------------------------------------
 do idij=1,nsploop

!  ------------------------------------------------------------------------
!  ----- Computation of Dij_so
!  ------------------------------------------------------------------------
   klmn1=1
   dijso(:,idij)=zero
   if (mod(idij,2)==1) then
     ispden=(1+idij)/2
     do klmn=1,lmn2_size
       if (indklmn(3,klmn)==0) then   ! il==jl
         klm=indklmn(1,klmn);kln=indklmn(2,klmn)
         ilm=indklmn(5,klmn);jlm=indklmn(6,klmn)
         fact=dijso_rad(kln);if (ilm>jlm) fact=-fact
         dijso(klmn1  ,idij)=fact*pawang%ls_ylm(1,klm,ispden)
         dijso(klmn1+1,idij)=fact*pawang%ls_ylm(2,klm,ispden)
       end if
       klmn1=klmn1+cplex_dij
     end do
   else if (idij==2) then
     do klmn=1,lmn2_size
       if (indklmn(3,klmn)==0) then   ! il==jl
         dijso(klmn1:klmn1+1,2)=-dijso(klmn1:klmn1+1,1)
       end if
       klmn1=klmn1+cplex_dij
     end do
   else if (idij==4) then
     do klmn=1,lmn2_size
       if (indklmn(3,klmn)==0) then   ! il==jl
         dijso(klmn1  ,4)=-dijso(klmn1  ,3)
         dijso(klmn1+1,4)= dijso(klmn1+1,3)
       end if
       klmn1=klmn1+cplex_dij
     end do
   end if

!  ----- End loop over idij
 end do

 LIBPAW_DEALLOCATE(dijso_rad)

end subroutine pawdijso
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawdiju
!! NAME
!! pawdiju
!!
!! FUNCTION
!! Compute the LDA+U contribution to the PAW pseudopotential strength Dij,
!! (for one atom only):
!!   Dijpawu^{\sigma}_{mi,ni,mj,nj}=
!!     \sum_{m,m'} [vpawu^{\sigma}_{m,m'}*phiphjint_{ni,nj}^{m,m'}]=
!!     [vpawu^{\sigma}_{mi,mj}*phiphjint_{ni,nj}]
!!
!! INPUTS
!!  cplex_dij=2 if dij is COMPLEX (as in the spin-orbit case), 1 if dij is REAL
!!  qphase=2 if dij contains a exp(-i.q.r) phase (as in the q<>0 RF case), 1 if not
!!  ndij= number of spin components
!!  nsppol=number of independent spin WF components
!!  pawtab <type(pawtab_type)>=paw tabulated starting data, for current atom
!!  vpawu(cplex_dij,lpawu*2+1,lpawu*2+1,ndij)=moments of LDA+U potential for current atom
!!  --- Optional arguments ---
!!    atvshift(natvshift,nsppol)=potential energy shift for lm channel & spin (current atom)
!!    fatvshift=factor that multiplies atvshift
!!    natvshift=number of atomic potential energy shifts (per atom)
!!
!! OUTPUT
!!  dijpawu(cplex_dij*qphase*lmn2_size,ndij)=  D_ij^U terms
!!    When Dij is complex (cplex_dij=2):
!!      dij(2*i-1,:) contains the real part, dij(2*i,:) contains the imaginary part
!!    When a exp(-i.q.r) phase is included (qphase=2):
!!      dij(1:cplex_dij*lmn2_size,:)
!!          contains the real part of the phase, i.e. D_ij*cos(q.r)
!!      dij(cplex_dij*lmn2_size+1:2*cplex_dij*lmn2_size,:)
!!          contains the imaginary part of the phase, i.e. D_ij*sin(q.r)
!!
!! PARENTS
!!      m_pawdij
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdiju(dijpawu,cplex_dij,qphase,ndij,nsppol,pawtab,vpawu,&
&                  natvshift,atvshift,fatvshift) ! optional arguments

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex_dij,ndij,nsppol,qphase
 integer,intent(in),optional :: natvshift
 real(dp),intent(in),optional :: fatvshift
!arrays
 real(dp),intent(out) :: dijpawu(:,:)
 real(dp),intent(in) :: vpawu(:,:,:,:)
 real(dp),intent(in),optional :: atvshift(:,:)
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: icount,idij,idijeff,idijend,im1,im2,in1,in2,klmn,klmn1,lmax,lmin,lmn2_size
 integer :: lpawu,natvshift_,nsploop
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: coeffpawu(:),dijpawu_idij(:),dijsymU(:,:)

! *************************************************************************

!Useful data
 lpawu=pawtab%lpawu
 lmn2_size=pawtab%lmn2_size
 natvshift_=0;if (present(natvshift)) natvshift_=natvshift

!Check data consistency
 if (qphase/=1) then
   msg='qphase=2 not available in pawdiju!'
   MSG_BUG(msg)
 end if
 if (size(dijpawu,1)/=cplex_dij*qphase*lmn2_size.or.size(dijpawu,2)/=ndij) then
   msg='invalid sizes for dijpawu !'
   MSG_BUG(msg)
 end if
 if (size(vpawu,1)/=cplex_dij.or.size(vpawu,2)/=2*lpawu+1.or.&
&    size(vpawu,3)/=2*lpawu+1.or.size(vpawu,4)/=ndij) then
   msg='invalid sizes for vpawu !'
   MSG_BUG(msg)
 end if
 if (natvshift_>0) then
   if ((.not.present(atvshift)).or.(.not.present(fatvshift))) then
     msg='when natvshift>0, atvshift and fatvshift arguments must be present !'
     MSG_BUG(msg)
   end if
   if (size(atvshift,1)/=natvshift.or.size(atvshift,2)/=nsppol) then
     msg='invalid sizes for atvshift !'
     MSG_BUG(msg)
   end if
 end if

!Init memory
 dijpawu=zero
 LIBPAW_ALLOCATE(dijpawu_idij,(cplex_dij*lmn2_size))
 LIBPAW_ALLOCATE(coeffpawu,(cplex_dij))
 if (ndij==4) then
   LIBPAW_ALLOCATE(dijsymU,(cplex_dij*lmn2_size,4))
 end if

!Loop over spin components
!----------------------------------------------------------
 nsploop=nsppol;if (ndij==4) nsploop=4
 do idij=1,nsploop
   if (idij<=nsppol.or.(ndij==4.and.idij<=3)) then

     idijend=idij+idij/3
     do idijeff=idij,idijend ! if ndij==4, idijeff is used to compute updn and dnup contributions

       dijpawu_idij=zero

!      Loop over (l,m,n) moments
!      ----------------------------------------------------------
       klmn1=1
       do klmn=1,lmn2_size
         im1=pawtab%klmntomn(1,klmn)
         im2=pawtab%klmntomn(2,klmn)
         lmin=pawtab%indklmn(3,klmn)
         lmax=pawtab%indklmn(4,klmn)

!        Select l=lpawu
         if (lmin==0.and.lmax==2*lpawu) then

!          Check consistency
           in1=pawtab%klmntomn(3,klmn)
           in2=pawtab%klmntomn(4,klmn)
           icount=in1+(in2*(in2-1))/2
           if (pawtab%ij_proj<icount)  then
             msg='LDA+U: Problem while computing dijexxc !'
             MSG_BUG(msg)
           end if

!          coeffpawu(:)=vpawu(:,im1,im2,idijeff) ! use real and imaginary part
           coeffpawu(:)=vpawu(:,im2,im1,idijeff) ! because of transposition in setnoccmmp (for the cplex_dij==2)

           if (natvshift_/=0.and.idij<3.and.im1==im2) then
             coeffpawu(1)=coeffpawu(1)+fatvshift*atvshift(im1,idij)
           end if
           if (cplex_dij==1) then   !cplex_dij=nspinor=1
             dijpawu_idij(klmn1)=pawtab%phiphjint(icount)*coeffpawu(1)
           elseif (cplex_dij==2) then   !cplex_dij=nspinor=2
             dijpawu_idij(klmn1  )=pawtab%phiphjint(icount)*coeffpawu(1)
             dijpawu_idij(klmn1+1)=pawtab%phiphjint(icount)*coeffpawu(2) !  spinor==2
           end if

         end if ! l selection
         klmn1=klmn1+cplex_dij
       end do ! klmn

       dijpawu(:,idij)=dijpawu_idij(:)
       if (ndij==4) dijsymU(:,idijeff)=dijpawu_idij(:)

     end do ! idijeff

   end if ! idij

   if (ndij==4.or.cplex_dij==2) then
     if (idij<=2)  then
       dijpawu(:,idij)=dijpawu(:,idij)
     else
       dijpawu(:,idij)=dijsymU(:,idij)
     end if
   end if

!End loop over spin components
!----------------------------------------------------------
 end do

!Free temporary memory spaces
 LIBPAW_DEALLOCATE(dijpawu_idij)
 LIBPAW_DEALLOCATE(coeffpawu)
 if (ndij==4) then
   LIBPAW_DEALLOCATE(dijsymU)
 end if

end subroutine pawdiju
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawdiju_euijkl
!! NAME
!! pawdiju_euijkl
!!
!! FUNCTION
!! Compute the LDA+U contribution to the PAW pseudopotential strength Dij (for one atom only).
!! Alternative to pawdiju using the following property:
!!     D_ij^pawu^{\sigma}_{mi,ni,mj,nj}=\sum_{k,l} [rho^{\sigma}_kl*e^U_ijkl]
!! The routine structure is very similar to the one of pawdijhartree.
!!
!! INPUTS
!!  cplex_dij=2 if dij is COMPLEX (as in the spin-orbit case), 1 if dij is REAL
!!  qphase=2 if dij contains a exp(-i.q.r) phase (as in the q<>0 RF case), 1 if not
!!  nspden=number of spin density components
!!  pawrhoij <type(pawrhoij_type)>= paw rhoij occupancies (and related data) for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data, for current atom
!!
!! OUTPUT
!!  diju(cplex_dij*qphase*lmn2_size,ndij)=  D_ij^U terms
!!  diju_im(cplex_dij*qphase*lmn2_size,ndij)= (see below)
!!    When Dij is complex (cplex_dij=2):
!!      dij(2*i-1,:) contains the real part, dij(2*i,:) contains the imaginary part
!!    When a exp(-i.q.r) phase is included (qphase=2):
!!      dij(1:cplex_dij*lmn2_size,:)
!!          contains the real part of the phase, i.e. D_ij*cos(q.r)
!!      dij(cplex_dij*lmn2_size+1:2*cplex_dij*lmn2_size,:)
!!          contains the imaginary part of the phase, i.e. D_ij*sin(q.r)
!!
!! NOTES
!! There are some subtleties :
!!   Contrary to eijkl, eu_ijkl is not invariant with respect to the permutation of i <--> j or k <--> l.
!!   So the correct expression of Dij is:
!!
!!     D_kl = sum_i<=j ( rho_ij eu_ijkl + (1-delta_ij) rho_ji eu_jikl )
!!
!!   In the following, we will use that: (according to the rules in pawpuxinit.F90)
!!    (a) eu_ijkl + eu_jikl =   eu_ijlk + eu_jilk (invariant      when exchanging k <--> l)
!!    (b) eu_ijkl - eu_jikl = - eu_ijlk + eu_jilk (change of sign when exchanging k <--> l)
!!   and :
!!    (c) eu_iikl = eu_iilk (if i=j, invariant when exchanging k <--> l)
!!    (d) eu_ijkk = eu_jikk (if k=l, invariant when exchanging i <--> j)
!!
!!   1) If qphase=1 (ipert=0 or q=0) we have simply:
!!        rho_ji = rho_ij^*
!!      So:
!!           D_kl  = sum_i<=j ( rho_ij eu_ijkl + (1-delta_ij) rho_ij^* eu_jikl )
!!      As eu_ijkl is real:
!! [I] :  Re(D_kl) = sum_i<=j Re(rho_ij) ( eu_ijkl + (1-delta_ij) eu_jikl )
!!        Im(D_kl) = sum_i<=j Im(rho_ij) ( eu_ijkl - (1-delta_ij) eu_jikl )
!!      So:
!!        Re(D_kl) = sum_i<=j Re(rho_ij) ( eu_ijlk + (1-delta_ij) eu_jilk ) =  Re(D_lk)  ( using (a) and (c) )
!!        Im(D_kl) = sum_i<=j Im(rho_ij) ( eu_ijlk - (1-delta_ij) eu_jilk ) = -Im(D_lk)  ( using (b) and (c) )
!!
!!   2) If qphase=2 (so ipert>0 and q/=0), we have:
!!        rho_ji = rhoA_ji + rhoB_ji
!!      where:
!!        rhoA_ji = rhoA_ij^*
!!        rhoB_ji = rhoB_ij
!!      So:
!!           D_kl = sum_i<=j ( rho_ij eu_ijkl + (1-delta_ij) (rhoA_ij^* + rhoB_ij) eu_jikl )
!!      As eu_ijkl is real:
!! [Ib] : Re(D_kl) = sum_i<=j Re(rho_ij)  ( eu_ijkl + (1-delta_ij) eu_jikl )  (same as [I])
!! [II] : Im(D_kl) = sum_i<=j Im(rhoB_ij) ( eu_ijkl + (1-delta_ij) eu_jikl )
!!                 + sum_i<=j Im(rhoA_ij) ( eu_ijkl - (1-delta_ij) eu_jikl )
!!      We note:
!!        Im(D_kl^A) = sum_i<=j Im(rhoA_ij) ( eu_ijkl - (1-delta_ij) eu_jikl )
!!        Im(D_kl^B) = sum_i<=j Im(rhoB_ij) ( eu_ijkl + (1-delta_ij) eu_jikl )
!!      We still have:
!!        Re(D_kl)  =  Re(D_lk)
!!      but:
!!        Im(D_kl^A) = -Im(D_lk^A)  ( using (b) and (c) )
!!        Im(D_kl^B) =  Im(D_lk^B)  ( using (a) and (c) )
!!
!! PARENTS
!!      m_pawdij,pawdenpot,pawdfptenergy
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdiju_euijkl(diju,cplex_dij,qphase,ndij,pawrhoij,pawtab,diju_im)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex_dij,ndij,qphase
!arrays
 real(dp),intent(out) :: diju(:,:)
 real(dp),intent(out),optional :: diju_im(:,:)
 type(pawrhoij_type),intent(in) :: pawrhoij
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: cplex_rhoij,iq,iq0_dij,iq0_rhoij,ilmn,ilmnp,irhoij,j0lmnp,jlmn,jlmnp,jrhoij
 integer :: klmn,klmnp,klmn1,lmn2_size,sig1,sig2
 logical :: compute_im
 character(len=500) :: msg
!arrays
 real(dp) :: ro(2)

! *************************************************************************

!Check data consistency
 lmn2_size=pawrhoij%lmn2_size
 if (size(diju,1)/=qphase*cplex_dij*lmn2_size.or.size(diju,2)/=ndij) then
   msg='invalid sizes for diju!'
   MSG_BUG(msg)
 end if
 if (pawrhoij%qphase<qphase) then
   msg='pawrhoij%qphase must be >=qphase!'
   MSG_BUG(msg)
 end if
 if (ndij==4) then
   msg='pawdiju_euijkl not yet available for ndij=4!'
   MSG_BUG(msg)
 end if

!Initialization
 diju=zero
 cplex_rhoij=pawrhoij%cplex_rhoij
 compute_im=(cplex_dij==2.and.cplex_rhoij==2)

!Loops over spin-components
 do sig2=1,max(pawrhoij%nspden,2)
   do sig1=1,max(pawrhoij%nspden,2)

     !Loop over phase exp(iqr) phase real/imaginary part
     do iq=1,qphase
       !First loop: we store the real part in dij(1 -> lmn2_size)
       !2nd loop: we store the imaginary part in dij(lmn2_size+1 -> 2*lmn2_size)
       iq0_dij=merge(0,cplex_dij*lmn2_size,iq==1)
       iq0_rhoij=merge(0,cplex_rhoij*lmn2_size,iq==1)

       !Loop over rhoij elements
       jrhoij=iq0_rhoij+1
       do irhoij=1,pawrhoij%nrhoijsel
         klmn=pawrhoij%rhoijselect(irhoij)
         ilmn=pawtab%indklmn(7,klmn)
         jlmn=pawtab%indklmn(8,klmn)

         ro(1:cplex_rhoij)=pawrhoij%rhoijp(jrhoij:jrhoij+cplex_rhoij-1,sig2)

         do jlmnp=1,pawtab%lmn_size
           j0lmnp=jlmnp*(jlmnp-1)/2
           do ilmnp=1,jlmnp
             klmnp=j0lmnp+ilmnp
             klmn1=iq0_dij+cplex_dij*(klmnp-1)+1

!            Re(D_kl) = sum_i<=j Re(rho_ij) ( eu_ijlk + (1-delta_ij) eu_jilk ) =  Re(D_lk)
             diju(klmn1,sig1)=diju(klmn1,sig1)+ro(1)*pawtab%euijkl(sig1,sig2,ilmn,jlmn,ilmnp,jlmnp)
             if (ilmn/=jlmn) diju(klmn1,sig1)=diju(klmn1,sig1)+ro(1)*pawtab%euijkl(sig1,sig2,jlmn,ilmn,ilmnp,jlmnp)

!            Im(D_kl) = sum_i<=j Im(rho_ij) ( eu_ijlk - (1-delta_ij) eu_jilk ) = -Im(D_lk)
             if (compute_im) then
               diju(klmn1+1,sig1)=diju(klmn1+1,sig1)+ro(2)*pawtab%euijkl(sig1,sig2,ilmn,jlmn,ilmnp,jlmnp)
               if (ilmn/=jlmn) diju(klmn1+1,sig1)=diju(klmn1+1,sig1)-ro(2)*pawtab%euijkl(sig1,sig2,jlmn,ilmn,ilmnp,jlmnp)
             end if

           end do ! k,l
         end do ! i,j

         jrhoij=jrhoij+cplex_rhoij
       end do

     end do ! q phase

   end do  !sig1
 end do !sig2


! OLD VERSION FROM LUCAS BAGUET
! ----------------------------------------------------------------
!  if (compute_diju_im) then
!    if (size(diju_im,1)/=lmn2_size.or.size(diju_im,2)/=ndij) then
!      msg='invalid sizes for diju_im !'
!      MSG_BUG(msg)
!    end if
!  end if
!  if (cplex_dij/=1) then
!    msg='pawdiju_euijkl not yet available for cplex_dij=2!'
!    MSG_ERROR(msg)
!  end if
!
!  lmn2_size=pawtab%lmn2_size
!  cplex_rhoij=pawrhoij%cplex_rhoij
!  compute_diju_im=(qphase==2.and.present(diju_im))
!
!  diju=zero
!  if (compute_diju_im) diju_im = zero
!
! !Real on-site quantities
!  if (qphase==1) then
!    do sig1=1,ndij
!      do sig2=1,ndij
!        jrhoij=1
!        do irhoij=1,pawrhoij%nrhoijsel
!          klmn=pawrhoij%rhoijselect(irhoij)
!          ilmn=pawtab%indklmn(7,klmn)
!          jlmn=pawtab%indklmn(8,klmn)
!          ro(1)=pawrhoij%rhoijp(jrhoij,sig2)
!          do jlmnp=1,pawtab%lmn_size
!            do ilmnp=1,jlmnp
!              klmn1 = ilmnp + jlmnp*(jlmnp-1)/2
!
! !            Thanks to Eq.[I] in the comment above:
!              diju(klmn1,sig1)=diju(klmn1,sig1)+ro(1)*pawtab%euijkl(sig1,sig2,ilmn,jlmn,ilmnp,jlmnp)
!              if (ilmn/=jlmn) then
!                diju(klmn1,sig1)=diju(klmn1,sig1)+ro(1)*pawtab%euijkl(sig1,sig2,jlmn,ilmn,ilmnp,jlmnp)
!              end if
!
!            end do
!          end do
!          jrhoij=jrhoij+cplex_rhoij
!        end do
!      end do
!    end do
!
! !Complex on-site quantities
!  else
!    do sig1=1,ndij
!      do sig2=1,ndij
!        jrhoij=1
!        do irhoij=1,pawrhoij%nrhoijsel
!          klmn=pawrhoij%rhoijselect(irhoij)
!          ilmn=pawtab%indklmn(7,klmn)
!          jlmn=pawtab%indklmn(8,klmn)
!          ro(1:2)=pawrhoij%rhoijp(jrhoij:jrhoij+1,sig2)
!          do jlmnp=1,pawtab%lmn_size
!            do ilmnp=1,jlmnp
!              klmn1 = ilmnp + jlmnp*(jlmnp-1)/2
!              kklmn1 = klmn1 + lmn2_size
!              ro_im = pawrhoij%rhoijim(klmn1,sig2)
!
! !            Thanks to Eq.[I] in the comment above:
!              diju(klmn1 ,sig1)=diju(klmn1 ,sig1)+ro(1)*pawtab%euijkl(sig1,sig2,ilmn,jlmn,ilmnp,jlmnp)
!              diju(kklmn1,sig1)=diju(kklmn1,sig1)+ro(2)*pawtab%euijkl(sig1,sig2,ilmn,jlmn,ilmnp,jlmnp)
!              diju(kklmn1,sig1)=diju(kklmn1,sig1)+ro_im*pawtab%euijkl(sig1,sig2,ilmn,jlmn,ilmnp,jlmnp)
!              if (compute_diju_im) then
!                diju_im(klmn1,sig1)=diju_im(klmn1,sig1)+ro_im*pawtab%euijkl(sig1,sig2,ilmn,jlmn,ilmnp,jlmnp)
!              end if
!
!              if (ilmn/=jlmn) then
!                diju(klmn1 ,sig1)=diju(klmn1 ,sig1)+ro(1)*pawtab%euijkl(sig1,sig2,jlmn,ilmn,ilmnp,jlmnp)
!                diju(kklmn1,sig1)=diju(klmn1 ,sig1)+ro(2)*pawtab%euijkl(sig1,sig2,jlmn,ilmn,ilmnp,jlmnp)
!                diju(kklmn1,sig1)=diju(kklmn1,sig1)-ro_im*pawtab%euijkl(sig1,sig2,jlmn,ilmn,ilmnp,jlmnp)
!                if (compute_diju_im) then
!                    diju_im(klmn1,sig1)=diju_im(klmn1,sig1)-ro_im*pawtab%euijkl(sig1,sig2,jlmn,ilmn,ilmnp,jlmnp)
!                end if
!              end if
!            end do
!          end do
!          jrhoij=jrhoij+cplex_rhoij
!        end do
!      end do
!    end do
!  end if

end subroutine pawdiju_euijkl
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawdijexxc
!! NAME
!! pawdijexxc
!!
!! FUNCTION
!! Compute the local Exact-Exchange contribution to the PAW pseudopotential strength Dij,
!! using a potential expressed as (l,m) spherical moments
!! (for one atom only; only for correlated electrons):
!!   D_ij^EXXC= < Phi_i|alpha*(VFock(correlated)-Vxc(n1_correlated)|Phi_j>
!!
!! INPUTS
!!  cplex_dij=2 if dij is COMPLEX (as in the spin-orbit case), 1 if dij is REAL
!!  qphase=2 if dij contains a exp(-i.q.r) phase (as in the q<>0 RF case), 1 if not
!!  lmselect(lm_size)=select the non-zero LM-moments of on-site potentials
!!  ndij= number of spin components
!!  nsppol=number of independent spin WF components
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data, for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data, for current atom
!!  vpawx(1,lmn2_size,ndij)=moments of exact exchange potential
!!                    for current atom and for correlated electrons
!!  vxc_ex(qphase*mesh_size,lm_size,nspden)=all-electron on-site XC potential for current atom
!!                    taken into account only valence correlated electrons
!!
!! OUTPUT
!!  dijexxc(cplex_dij*lmn2_size,ndij)=  D_ij^Exact-Exchange terms
!!    When Dij is complex (cplex_dij=2):
!!      dij(2*i-1,:) contains the real part, dij(2*i,:) contains the imaginary part
!!    When a exp(-i.q.r) phase is included (qphase=2):
!!      dij(1:cplex_dij*lmn2_size,:)
!!          contains the real part of the phase, i.e. D_ij*cos(q.r)
!!      dij(cplex_dij*lmn2_size+1:2*cplex_dij*lmn2_size,:)
!!          contains the imaginary part of the phase, i.e. D_ij*sin(q.r)
!!
!! PARENTS
!!      m_pawdij
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijexxc(dijexxc,cplex_dij,qphase,lmselect,ndij,nspden,nsppol,&
&                     pawang,pawrad,pawtab,vpawx,vxc_ex)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex_dij,ndij,nspden,nsppol,qphase
 type(pawang_type),intent(in) :: pawang
!arrays
 logical :: lmselect(:)
 real(dp),intent(in) :: vpawx(:,:,:),vxc_ex(:,:,:)
 real(dp),intent(out) :: dijexxc(:,:)
 type(pawrad_type),intent(in) :: pawrad
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: icount,idij,idijend,ij_size,iln,in1,in2,ir,ir1,isel,ispden,ivxc
 integer :: jln,j0ln,klm,klm1,klmn,klmn1,klmn2,kln,lexexch,ln_min,ln_max,lmax,lmin
 integer :: lm_size,lmn2_size,mesh_size,nsploop
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: dijexxc_idij(:),ff(:),gg(:),vxcij1(:)

! *************************************************************************

!Useful data
 lm_size=pawtab%lcut_size**2
 lmn2_size=pawtab%lmn2_size
 ij_size=pawtab%ij_size
 mesh_size=pawtab%mesh_size
 lexexch=pawtab%lexexch
 ln_min=pawtab%lnproju(1)
 ln_max=pawtab%lnproju(pawtab%nproju)

!Check data consistency
 if (qphase==2) then
   msg='pawdijexx not available for qphase=2!'
   MSG_BUG(msg)
 end if
 if (size(dijexxc,1)/=cplex_dij*qphase*lmn2_size.or.size(dijexxc,2)/=ndij) then
   msg='invalid sizes for dijexxc!'
   MSG_BUG(msg)
 end if
 if (size(lmselect)/=lm_size) then
   msg='invalid size for lmselect!'
   MSG_BUG(msg)
 end if
 if (size(vxc_ex,1)/=qphase*mesh_size.or.size(vxc_ex,2)/=lm_size.or.&
&    size(vxc_ex,3)/=nspden) then
   msg='invalid sizes for vxc_ex!'
   MSG_BUG(msg)
 end if
 if (size(vpawx,1)/=1.or.size(vpawx,2)/=lmn2_size.or.&
&    size(vpawx,3)/=ndij) then
   msg='invalid sizes for vpawx!'
   MSG_BUG(msg)
 end if

!Init memory
 dijexxc=zero
 LIBPAW_ALLOCATE(dijexxc_idij,(qphase*lmn2_size))
 LIBPAW_ALLOCATE(vxcij1,(qphase*ij_size))
 LIBPAW_ALLOCATE(ff,(mesh_size))
 LIBPAW_ALLOCATE(gg,(mesh_size))

!----------------------------------------------------------
!Loop over spin components
!----------------------------------------------------------
 nsploop=nsppol;if (ndij==4) nsploop=4
 do idij=1,nsploop

   if (idij<=nsppol.or.(ndij==4.and.idij<=3)) then

     idijend=idij+idij/3
     do ispden=idij,idijend

       dijexxc_idij=zero

       ivxc=ispden
       !Take into account nspden=1/nspinor=2 case
       if (ndij/=nspden.and.ispden==2) ivxc=1
       if (ndij/=nspden.and.ispden> 2) cycle

!      ----------------------------------------------------------
!      Summing over (l,m) moments
!      ----------------------------------------------------------
       do klm=1,lm_size
         if (lmselect(klm)) then

!          ===== Vxc_ij_1 (tmp) =====
           vxcij1=zero
           if (qphase==1) then
             do jln=ln_min,ln_max
               j0ln=jln*(jln-1)/2
               do iln=ln_min,jln
                 kln=j0ln+iln
                 ff(1:mesh_size)= &
&                  vxc_ex(1:mesh_size,klm,ivxc)*pawtab%phiphj(1:mesh_size,kln)
                 call simp_gen(vxcij1(kln),ff,pawrad)
               end do
             end do
           else
             do jln=ln_min,ln_max
               j0ln=jln*(jln-1)/2
               do iln=ln_min,jln
                 kln=j0ln+iln
                 do ir=1,mesh_size
                   ir1=2*ir
                   ff(ir)= &
&                   vxc_ex(ir1-1,klm,ivxc)*pawtab%phiphj(ir,kln)
                   gg(ir)= &
&                   vxc_ex(ir1,klm,ivxc)*pawtab%phiphj(ir,kln)
                 end do
                 call simp_gen(vxcij1(2*kln-1),ff,pawrad)
                 call simp_gen(vxcij1(2*kln  ),gg,pawrad)
               end do
             end do
           end if

!          ===== Accumulate Vxc_ij_1 over klm moments =====
           if (qphase==1) then
             do klmn=1,lmn2_size
               lmin=pawtab%indklmn(3,klmn)
               lmax=pawtab%indklmn(4,klmn)
               if (lmin==0.and.lmax==2*lexexch) then
                 klm1=pawtab%indklmn(1,klmn)
                 kln=pawtab%indklmn(2,klmn)
                 isel=pawang%gntselect(klm,klm1)
                 if (isel>0) dijexxc_idij(klmn)=dijexxc_idij(klmn) &
&                                  +vxcij1(kln)*pawang%realgnt(isel)
               end if
             end do ! Loop klmn
           else ! qphase==2
             klmn1=1
             do klmn=1,lmn2_size
               lmin=pawtab%indklmn(3,klmn)
               lmax=pawtab%indklmn(4,klmn)
               if (lmin==0.and.lmax==2*lexexch) then
                 klm1=pawtab%indklmn(1,klmn)
                 kln=pawtab%indklmn(2,klmn)
                 isel=pawang%gntselect(klm,klm1)
                 if (isel>0) then
                   dijexxc_idij(klmn1  )=dijexxc_idij(klmn1) &
&                                     +vxcij1(2*kln-1)*pawang%realgnt(isel)
                   dijexxc_idij(klmn1+1)=dijexxc_idij(klmn1+1) &
&                                     +vxcij1(2*kln  )*pawang%realgnt(isel)
                 end if
               end if
               klmn1=klmn1+qphase
             end do ! Loop klmn
           end if

         end if ! lmselect
       end do  ! Loop klm

!      Mix Hartree and GGA terms
       if (qphase==1) then
         do klmn=1,lmn2_size
           lmin=pawtab%indklmn(3,klmn)
           lmax=pawtab%indklmn(4,klmn)
           if (lmin==0.and.lmax==2*lexexch) then
             in1=pawtab%klmntomn(3,klmn)
             in2=pawtab%klmntomn(4,klmn)
             icount=in1+(in2*(in2-1))/2
             if(pawtab%ij_proj<icount)  then
               msg='PAW local exact-exchange: Problem while computing dijexxc !'
               MSG_BUG(msg)
             end if
             dijexxc_idij(klmn)=pawtab%exchmix &
&                              *(vpawx(1,klmn,idij)-dijexxc_idij(klmn))
           end if
         end do
       else ! qphase=2
         klmn1=1
         do klmn=1,lmn2_size
           lmin=pawtab%indklmn(3,klmn)
           lmax=pawtab%indklmn(4,klmn)
           if (lmin==0.and.lmax==2*lexexch) then
             in1=pawtab%klmntomn(3,klmn)
             in2=pawtab%klmntomn(4,klmn)
             icount=in1+(in2*(in2-1))/2
             if(pawtab%ij_proj<icount)  then
               msg='PAW local exact-exchange: Problem while computing dijexxc !'
               MSG_BUG(msg)
             end if
             dijexxc_idij(klmn1)  =pawtab%exchmix &
&                                 *(vpawx(1,klmn,idij)-dijexxc_idij(klmn1))
             dijexxc_idij(klmn1+1)=pawtab%exchmix &
&                                 *(vpawx(1,klmn,idij)-dijexxc_idij(klmn1+1))
           end if
           klmn1=klmn1+qphase
         end do ! Loop klmn
       end if

!      ----------------------------------------------------------
!      Deduce some part of Dij according to symmetries
!      ----------------------------------------------------------

       !if ispden=1 => real part of D^11_ij
       !if ispden=2 => real part of D^22_ij
       !if ispden=3 => real part of D^12_ij
       !if ispden=4 => imaginary part of D^12_ij
       klmn1=max(1,ispden-2);klmn2=1
       do klmn=1,lmn2_size
         dijexxc(klmn1,idij)=dijexxc_idij(klmn2)
         klmn1=klmn1+cplex_dij
         klmn2=klmn2+qphase
       end do
       if (qphase==2) then
         !Same storage with exp^(-i.q.r) phase
         klmn1=max(1,ispden-2)+lmn2_size*cplex_dij;klmn2=2
         do klmn=1,lmn2_size
           dijexxc(klmn1,idij)=dijexxc_idij(klmn2)
           klmn1=klmn1+cplex_dij
           klmn2=klmn2+qphase
         end do
       endif

     end do !ispden

   !Non-collinear: D_ij(:,4)=D^21_ij=D^12_ij^*
   else if (nspden==4.and.idij==4) then
     dijexxc(:,idij)=dijexxc(:,idij-1)
     if (cplex_dij==2) then
       do klmn=2,lmn2_size*cplex_dij,cplex_dij
         dijexxc(klmn,idij)=-dijexxc(klmn,idij)
       end do
       if (qphase==2) then
         do klmn=2+lmn2_size*cplex_dij,2*lmn2_size*cplex_dij,cplex_dij
           dijexxc(klmn,idij)=-dijexxc(klmn,idij)
         end do
       end if
     end if

   !Antiferro: D_ij(:,2)=D^down_ij=D^up_ij
   else if (nsppol==1.and.idij==2) then
     dijexxc(:,idij)=dijexxc(:,idij-1)
   end if

!----------------------------------------------------------
!End loop on spin density components
 end do

!Free temporary memory spaces
 LIBPAW_DEALLOCATE(dijexxc_idij)
 LIBPAW_DEALLOCATE(vxcij1)
 LIBPAW_DEALLOCATE(ff)
 LIBPAW_DEALLOCATE(gg)

end subroutine pawdijexxc
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawdijfr
!!
!! NAME
!! pawdijfr
!!
!! FUNCTION
!! PAW, Response Function only:
!!      Compute frozen part of psp strength Dij due to 1st-order compensation density
!!      and first order local potential:
!!      Dijfr    =Int_R^3{vtrial*Sum_LM[Q_ij_q^LM^(1)] + Vloc^(1)*Sum_LM[Q_ij_q^LM]}
!!      Depends on q wave vector but not on first-order wave-function.
!!
!! INPUTS
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!  idir=direction of atomic displacement (in case of phonons perturb.)
!!  ipert=index of perturbation
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  mpi_comm_grid=--optional-- MPI communicator over real space grid components
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components
!!  nsppol=number of independent spin WF components
!!  ntypat=number of types of atoms
!!  option=0: computes full frozen part of Dij
!!         1: computes frozen part of Dij without contribution from Vpsp1
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(my_natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  qphase=2 if dij contains a exp(-i.q.r) phase (as in the q<>0 RF case), 1 if not
!!  qphon(3)=wavevector of the phonon
!!  rprimd(3,3)=dimensional primitive translations for real space
!!  ucvol=unit cell volume (bohr^3)
!!  vpsp1(qphase*nfft)= first-order change of local potential
!!  vtrial(nfft,nspden)= total GS potential
!!  vxc(nfft,nspden)=XC potential
!!  xred(3,my_natom)= reduced atomic coordinates
!!
!! OUTPUT
!!  paw_ij1(iatom)%dijfr(cplex_dij*qphase*lmn2_size,nspden)=
!!                  frozen contribution to psp strength Dij
!!                  =Int_R^3{vtrial*Sum_LM[Q_ij_q^LM^(1)] + Vloc^(1)*Sum_LM[Q_ij_q^LM]}
!!    When Dij is complex (cplex_dij=2):
!!      dij(2*i-1,:) contains the real part, dij(2*i,:) contains the imaginary part
!!    When a exp(-i.q.r) phase is included (qphase=2):
!!      dij(1:cplex_dij*lmn2_size,:)
!!          contains the real part of the phase, i.e. D_ij*cos(q.r)
!!      dij(cplex_dij*lmn2_size+1:2*cplex_dij*lmn2_size,:)
!!          contains the imaginary part of the phase, i.e. D_ij*sin(q.r)
!!
!! PARENTS
!!      d2frnl,dfpt_nstpaw,dfpt_rhofermi,dfpt_scfcv
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijfr(gprimd,idir,ipert,my_natom,natom,nfft,ngfft,nspden,nsppol,ntypat,&
&          option,paw_ij1,pawang,pawfgrtab,pawrad,pawtab,qphase,qphon,rprimd,ucvol,&
&          vpsp1,vtrial,vxc,xred,&
&          mpi_atmtab,comm_atom,mpi_comm_grid) ! optional arguments (parallelism)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert,my_natom,natom,nfft,nspden,nsppol,ntypat,option,qphase
 integer,optional,intent(in) :: comm_atom,mpi_comm_grid
 real(dp),intent(in) :: ucvol
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gprimd(3,3),qphon(3),rprimd(3,3)
 real(dp),intent(in) :: vpsp1(qphase*nfft),vtrial(nfft,nspden),vxc(nfft,nspden)
 real(dp),intent(in) :: xred(3,natom)
 type(paw_ij_type),intent(inout) :: paw_ij1(my_natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: cplex_dij,cplex_nspden,dplex_nsp,dplex_q,iatom,iatom_tot,ic,idij,idijend,ier,ils,ilslm,isel
 integer :: ispden,istr,itypat,jc,klm,klmn,klmn1,klmn2,kln,lm_size,lmn2_size,lm0,lmax,lmin,mesh_size
 integer :: mm,my_comm_atom,my_comm_grid,mu,mua,mub,ndij,nfftot,nfgd,nsploop
 integer :: optgr0,optgr1,optgr2,usexcnhat
 logical :: has_qphase,my_atmtab_allocated,need_dijfr_1,need_dijfr_2,need_dijfr_3,need_dijfr_4
 logical :: paral_atom,qne0,testdij1,testdij2,testdij3
 real(dp) :: c1,fact,intg,rg1
 character(len=500) :: msg
!arrays
 integer,parameter :: m_index(3)=(/1,-1,0/)
 integer,pointer :: my_atmtab(:)
 integer,parameter :: alpha(9)=(/1,2,3,3,3,2,2,1,1/),beta(9)=(/1,2,3,2,1,1,3,3,2/)
 real(dp) :: contrib(2)
 real(dp),allocatable :: ff(:),intv(:,:),intv1(:,:),intv2(:,:),intvloc(:,:),intv_tmp(:,:)
 real(dp),allocatable :: rg(:),vloc(:,:)

! *************************************************************************

!Nothing to be done for DDK
 if (ipert==natom+1.or.ipert==natom+10) return

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 my_comm_grid=xmpi_comm_self;if (present(mpi_comm_grid)) my_comm_grid=mpi_comm_grid
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Compatibility tests
 qne0=(qphon(1)**2+qphon(2)**2+qphon(3)**2>=1.d-15)
 if (my_natom>0) then
   if (qne0.and.qphase==1) then
     msg='qphase must be 2 when q<>0!'
     MSG_BUG(msg)
   end if
   if (paw_ij1(1)%qphase/=qphase) then
     msg='paw_ij1()%qphase and qphase must be equal !'
     MSG_BUG(msg)
   end if
   if (paw_ij1(1)%has_dijfr==0) then
     msg='pawdij1()%dijfr must be allocated !'
     MSG_BUG(msg)
   end if
   testdij1=(ipert<=natom.and.option==0.and.pawfgrtab(1)%gylm_allocated==0)
   testdij2=(ipert<=natom.and.pawfgrtab(1)%gylmgr_allocated==0)
   testdij3=(testdij2.and.qne0.and.pawfgrtab(1)%expiqr_allocated==0)
   if ((testdij1.or.testdij2.or.testdij3).and.pawfgrtab(1)%rfgd_allocated==0) then
     msg='pawfgrtab()%rfgd array must be allocated  !'
     MSG_BUG(msg)
   end if
 end if

!Get correct index of strain pertubation
 if (ipert==natom+3) istr = idir
 if (ipert==natom+4) istr = idir + 3

!Some inits
 usexcnhat=maxval(pawtab(1:ntypat)%usexcnhat)
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 fact=ucvol/dble(nfftot)
 cplex_nspden=merge(1,2,nspden/=4)
 dplex_nsp=cplex_nspden-1
 dplex_q=qphase-1

!Loops over  atoms
 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

!  Select which part of Dijfr to compute
   need_dijfr_1=(ipert==iatom_tot.and.paw_ij1(iatom)%has_dijfr==1)
   need_dijfr_2=(ipert<=natom.and.paw_ij1(iatom)%has_dijfr==1.and.(option==0))
   need_dijfr_3=((ipert==natom+2.or.ipert==natom+11).and.paw_ij1(iatom)%has_dijfr==1)
   need_dijfr_4=((ipert==natom+3.or.ipert==natom+4).and.paw_ij1(iatom)%has_dijfr==1)

   if ((.not.need_dijfr_1).and.(.not.need_dijfr_2).and.(.not.need_dijfr_3).and.(.not.need_dijfr_4)) then
     if (paw_ij1(iatom)%has_dijfr>0) then
       paw_ij1(iatom)%dijfr=zero ; paw_ij1(iatom)%has_dijfr=2
     end if
     cycle
   end if

!  Some atom-dependent quantities
   itypat=pawfgrtab(iatom)%itypat
   lm_size=pawtab(itypat)%lcut_size**2
   lmn2_size=pawtab(itypat)%lmn2_size
   cplex_dij=paw_ij1(iatom)%cplex_dij
   ndij=paw_ij1(iatom)%ndij

!  Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
   nfgd=0
   if (need_dijfr_1.or.need_dijfr_2.or.need_dijfr_4) then
     nfgd=pawfgrtab(iatom)%nfgd
     if (((need_dijfr_2.or.need_dijfr_4).and.(pawfgrtab(iatom)%gylm_allocated==0)).or.&
&     ((need_dijfr_1).and.(pawfgrtab(iatom)%gylmgr_allocated==0))) then
       optgr0=0;optgr1=0;optgr2=0
       if ((need_dijfr_2.or. need_dijfr_4).and.(pawfgrtab(iatom)%gylm_allocated==0)) then
         if (allocated(pawfgrtab(iatom)%gylm))  then
           LIBPAW_DEALLOCATE(pawfgrtab(iatom)%gylm)
         end if
         LIBPAW_ALLOCATE(pawfgrtab(iatom)%gylm,(nfgd,lm_size))
         pawfgrtab(iatom)%gylm_allocated=2;optgr0=1
       end if
       if ((need_dijfr_1.or.need_dijfr_4).and.(pawfgrtab(iatom)%gylmgr_allocated==0)) then
         if (allocated(pawfgrtab(iatom)%gylmgr))  then
           LIBPAW_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
         end if
         LIBPAW_ALLOCATE(pawfgrtab(iatom)%gylmgr,(3,nfgd,lm_size))
         pawfgrtab(iatom)%gylmgr_allocated=2;optgr1=1
       end if
       if (optgr0+optgr1+optgr2>0) then
         call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,pawfgrtab(iatom)%gylmgr2,&
&             lm_size,nfgd,optgr0,optgr1,optgr2,pawtab(itypat),pawfgrtab(iatom)%rfgd)
       end if
     end if
   end if

!  Eventually compute exp(-i.q.r) factors for the current atom (if not already done)
   has_qphase=(qne0.and.qphase==2)
   if (need_dijfr_2) then
     if (has_qphase.and.(pawfgrtab(iatom)%expiqr_allocated==0)) then
       if (allocated(pawfgrtab(iatom)%expiqr))  then
         LIBPAW_DEALLOCATE(pawfgrtab(iatom)%expiqr)
       end if
       LIBPAW_ALLOCATE(pawfgrtab(iatom)%expiqr,(2,nfgd))
       call pawexpiqr(pawfgrtab(iatom)%expiqr,gprimd,nfgd,qphon,&
&                     pawfgrtab(iatom)%rfgd,xred(:,iatom_tot))
       pawfgrtab(iatom)%expiqr_allocated=2
     end if
     has_qphase=(pawfgrtab(iatom)%expiqr_allocated/=0)
   end if

!  Loop over spin components
   nsploop=nsppol;if (ndij==4) nsploop=4
   do idij=1,nsploop
     if (idij<=nsppol.or.(nspden==4.and.idij<=3)) then

       idijend=idij+idij/3
       do ispden=idij,idijend

         LIBPAW_ALLOCATE(intv,(qphase*cplex_nspden,lm_size))
         intv(:,:) = zero

!        ============ Phonons ====================================
         if (ipert<=natom) then

           if (need_dijfr_1.or.need_dijfr_2) then

             LIBPAW_ALLOCATE(intv1,(cplex_nspden,lm_size))
             LIBPAW_ALLOCATE(intv2,(qphase,lm_size))
             intv1(:,:)=zero ; intv2(:,:)=zero

!            First part: Int_R^3{vtrial*Sum_LM[Q_ij_q^LM^(1)]}
             if (need_dijfr_1) then

!              ----- Retrieve potential Vlocal (subtle if nspden=4 ;-)
               LIBPAW_ALLOCATE(vloc,(cplex_nspden,nfgd))
               if (nspden/=4) then
                 if (usexcnhat==0) then
                   do ic=1,nfgd
                     jc=pawfgrtab(iatom)%ifftsph(ic)
                     vloc(1,ic)=vtrial(jc,ispden)-vxc(jc,ispden)
                   end do
                 else
                   do ic=1,nfgd
                     vloc(1,ic)=vtrial(pawfgrtab(iatom)%ifftsph(ic),ispden)
                   end do
                 end if
               else ! nspden==4
                 if (ispden<=2) then
                   if (usexcnhat==0) then
                     do ic=1,nfgd
                       jc=pawfgrtab(iatom)%ifftsph(ic)
                       vloc(1,ic)=vtrial(jc,ispden)-vxc(jc,ispden)
                       vloc(2,ic)=zero
                     end do
                   else
                     do ic=1,nfgd
                       jc=pawfgrtab(iatom)%ifftsph(ic)
                       vloc(1,ic)=vtrial(jc,ispden)
                       vloc(2,ic)=zero
                     end do
                   end if
                 else if (ispden==3) then
                   if (usexcnhat==0) then
                     vloc(:,:)=zero
                   else
                     do ic=1,nfgd
                       jc=pawfgrtab(iatom)%ifftsph(ic)
                       vloc(1,ic)=vtrial(jc,3)
                       vloc(2,ic)=vtrial(jc,4)
                     end do
                   end if
                 else ! ispden=4
                   vloc(2,1:nfgd)=-vloc(2,1:nfgd)
                 end if
               end if

!              ----- Compute Integral [ Vtrial(r).(g_l(r).Y_lm(r))^(1) dr ]
               LIBPAW_ALLOCATE(intv_tmp,(cplex_nspden,3))
               do ilslm=1,lm_size
                 intv_tmp=zero
                 do ic=1,nfgd
                   do mu=1,3
!                    Minus sign because dg(r-R)/dR = -dg(r-R)/dr
                     contrib(1:cplex_nspden)=-vloc(1:cplex_nspden,ic)*pawfgrtab(iatom)%gylmgr(mu,ic,ilslm)
                     intv_tmp(1:cplex_nspden,mu)=intv_tmp(1:cplex_nspden,mu)+contrib(1:cplex_nspden)
                   end do
                 end do
!                Convert from cartesian to reduced coordinates
                 intv1(1:cplex_nspden,ilslm)=intv1(1:cplex_nspden,ilslm) &
&                   +(rprimd(1,idir)*intv_tmp(1:cplex_nspden,1) &
&                    +rprimd(2,idir)*intv_tmp(1:cplex_nspden,2) &
&                    +rprimd(3,idir)*intv_tmp(1:cplex_nspden,3))
               end do
               LIBPAW_DEALLOCATE(vloc)
               LIBPAW_DEALLOCATE(intv_tmp)
             end if ! need_dijfr_1

!            2nd part: Int_R^3{Vloc^(1)*Sum_LM[Q_ij_q^LM]}
             if (need_dijfr_2) then

               if (ispden==1) then

!                ----- Retrieve potential Vloc^(1)
                 LIBPAW_ALLOCATE(vloc,(qphase,nfgd))
                 do ic=1,nfgd
                   jc=qphase*pawfgrtab(iatom)%ifftsph(ic)-dplex_q
                   vloc(1:qphase,ic)=vpsp1(jc:jc+dplex_q)
                 end do

!                ----- Compute Integral [ Vloc^(1)(r).g_l(r).Y_lm(r) ]
                 LIBPAW_ALLOCATE(intvloc,(qphase,lm_size))
                 intvloc=zero
                 if (has_qphase) then
                   if (qphase==1) then
                     do ilslm=1,lm_size
                       do ic=1,nfgd
                         contrib(1)=vloc(1,ic)*pawfgrtab(iatom)%gylm(ic,ilslm)
                         intvloc(1,ilslm)=intvloc(1,ilslm)+contrib(1)*pawfgrtab(iatom)%expiqr(1,ic)
                       end do
                     end do
                   else
                     do ilslm=1,lm_size
                       do ic=1,nfgd
                         contrib(1:2)=vloc(1:2,ic)*pawfgrtab(iatom)%gylm(ic,ilslm)
                         intvloc(1,ilslm)=intvloc(1,ilslm)+contrib(1)*pawfgrtab(iatom)%expiqr(1,ic) &
&                                                         -contrib(2)*pawfgrtab(iatom)%expiqr(2,ic)
                         intvloc(2,ilslm)=intvloc(2,ilslm)+contrib(1)*pawfgrtab(iatom)%expiqr(2,ic) &
&                                                         +contrib(2)*pawfgrtab(iatom)%expiqr(1,ic)
                       end do
                     end do
                   end if
                 else ! no phase
                   do ilslm=1,lm_size
                     do ic=1,nfgd
                       contrib(1:qphase)=vloc(1:qphase,ic)*pawfgrtab(iatom)%gylm(ic,ilslm)
                       intvloc(1:qphase,ilslm)=intvloc(1:qphase,ilslm)+contrib(1:qphase)
                     end do
                   end do
                 end if
                 LIBPAW_DEALLOCATE(vloc)
               end if ! ispden=1

               !Add to previous contribution
               if (ispden<=min(nspden,2)) then
                 intv2(1:qphase,1:lm_size)=intv2(1:qphase,1:lm_size)+intvloc(1:qphase,1:lm_size)
                 if (ispden==min(nspden,2)) then
                   LIBPAW_DEALLOCATE(intvloc)
                 end if
               end if
             end if ! need_dijfr_2

!            Sum contributions and apply ucvol/nfft factor on integral
             intv(1:cplex_nspden,1:lm_size)=intv1(1:cplex_nspden,1:lm_size)
             intv(1,1:lm_size)=intv(1,1:lm_size)+intv2(1,1:lm_size)
             if (qphase==2) intv(cplex_nspden+1,1:lm_size)=intv(cplex_nspden+1,1:lm_size)+intv2(2,1:lm_size)
             intv(:,:)=fact*intv(:,:)
             LIBPAW_DEALLOCATE(intv1)
             LIBPAW_DEALLOCATE(intv2)

!            --- Reduction in case of parallelization ---
             call xmpi_sum(intv,my_comm_grid,ier)

             paw_ij1(iatom)%dijfr(:,ispden)=zero

!            ---- Loop over (i,j) components
             klmn1=1;klmn2=1+lmn2_size*cplex_dij
             do klmn=1,lmn2_size
               klm =pawtab(itypat)%indklmn(1,klmn)
               lmin=pawtab(itypat)%indklmn(3,klmn)
               lmax=pawtab(itypat)%indklmn(4,klmn)
               do ils=lmin,lmax,2
                 lm0=ils**2+ils+1
                 do mm=-ils,ils
                   ilslm=lm0+mm;isel=pawang%gntselect(lm0+mm,klm)
                   if (isel>0) then
                     !The following works only because cplex_nspden<=cplex_dij
                     paw_ij1(iatom)%dijfr(klmn1:klmn1+dplex_nsp,ispden)= &
    &                 paw_ij1(iatom)%dijfr(klmn1:klmn1+dplex_nsp,ispden) &
    &                 +pawtab(itypat)%qijl(ilslm,klmn)*intv(1:cplex_nspden,ilslm)
                     if (qphase==2) then
                       paw_ij1(iatom)%dijfr(klmn2:klmn2+dplex_nsp,ispden)= &
    &                   paw_ij1(iatom)%dijfr(klmn2:klmn2+dplex_nsp,ispden) &
    &                   +pawtab(itypat)%qijl(ilslm,klmn)*intv(1+cplex_nspden:2*cplex_nspden,ilslm)
                     end if
                   end if
                 end do
               end do
               klmn1=klmn1+cplex_dij;klmn2=klmn2+cplex_dij
             end do

!            Dijfr is marked as computed
             paw_ij1(iatom)%has_dijfr=2

           end if

!        ============ Electric field perturbation =======================
         else if (ipert==natom+2.or.ipert==natom+11) then

           if (need_dijfr_3) then

!            The following factor arises in expanding the angular dependence of the dipole
!            vector in terms of real spherical harmonics. The real spherical harmonics are as
!            in the routine initylmr.F90;
!            see http://www.unioviedo.es/qcg/art/Theochem419-19-ov-BF97-rotation-matrices.pdf
             c1 = sqrt(four_pi/three)
             mesh_size=pawtab(itypat)%mesh_size

             if (ispden==1) then

               LIBPAW_ALLOCATE(ff,(mesh_size))
               LIBPAW_ALLOCATE(rg,(3))

!              loop over basis state pairs for this atom
               klmn1=1
               do klmn = 1, paw_ij1(iatom)%lmn2_size
                 klm =pawtab(itypat)%indklmn(1,klmn)
                 kln =pawtab(itypat)%indklmn(2,klmn)
                 lmin=pawtab(itypat)%indklmn(3,klmn)
                 lmax=pawtab(itypat)%indklmn(4,klmn)

!                Select only l=1, because the dipole is a vector operator
                 if (lmin==1) then
                   lm0=3  ! (l^2+l+1) for l=1

!                  Computation of <phi_i|r|phi_j>- <tphi_i|r|tphi_j>
!                  the dipole vector has radial dependence r
                   ff(1:mesh_size)=(pawtab(itypat)%phiphj(1:mesh_size,kln)&
&                   -pawtab(itypat)%tphitphj(1:mesh_size,kln))&
&                   *pawrad(itypat)%rad(1:mesh_size)
!                   call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
                   call simp_gen(intg,ff,pawrad(itypat))

!                  Compute <S_li_mi|r-R|S_lj_mj>: use a real Gaunt expression (with selection rule)
                   rg(1:3)=zero
                   do ic=1,3
                     isel=pawang%gntselect(lm0+m_index(ic),klm)
                     if (isel>0) rg(ic)=pawang%realgnt(isel)
                   end do

!                  Translate from cartesian to reduced coordinates (in idir direction)
                   rg1=gprimd(1,idir)*rg(1)+gprimd(2,idir)*rg(2)+gprimd(3,idir)*rg(3)

!                  Build sqrt(4pi/3).<S_li_mi|r-R|S_lj_mj>.(<phi_i|r-R|phi_j>- <tphi_i|r-R|tphi_j>
                   paw_ij1(iatom)%dijfr(klmn1,ispden)=c1*rg1*intg
                   if (cplex_dij==2) paw_ij1(iatom)%dijfr(klmn1+1,ispden)=zero

                 else
                   paw_ij1(iatom)%dijfr(klmn1,ispden)=zero
                 end if ! end gaunt constraint

                 klmn1=klmn1+cplex_dij
               end do ! end loop over lmn2_size pairs of basis states
               LIBPAW_DEALLOCATE(ff)
               LIBPAW_DEALLOCATE(rg)

!            Dijfr is spin-independent for electric field case
             else if (ispden==2) then
               paw_ij1(iatom)%dijfr(:,ispden)=paw_ij1(iatom)%dijfr(:,1)
             else
               paw_ij1(iatom)%dijfr(:,ispden)=zero
             end if

!            Dijfr is marked as computed
             paw_ij1(iatom)%has_dijfr=2
           end if

!        ============ Elastic tensor ===============================
         else if (ipert==natom+3.or.ipert==natom+4) then

!          ----- Retrieve potential Vlocal (subtle if nspden=4 ;-)
           LIBPAW_ALLOCATE(vloc,(cplex_nspden,nfgd))
           if (nspden/=4) then
             if (usexcnhat==0) then
               do ic=1,nfgd
                 jc=pawfgrtab(iatom)%ifftsph(ic)
                 vloc(1,ic)=vtrial(jc,ispden)-vxc(jc,ispden)
               end do
             else
               do ic=1,nfgd
                 vloc(1,ic)=vtrial(pawfgrtab(iatom)%ifftsph(ic),ispden)
               end do
             end if
           else ! nspden/=4
             if (ispden<=2) then
               if (usexcnhat==0) then
                 do ic=1,nfgd
                   jc=pawfgrtab(iatom)%ifftsph(ic)
                   vloc(1,ic)=vtrial(jc,ispden)-vxc(jc,ispden)
                   vloc(2,ic)=zero
                 end do
               else
                 do ic=1,nfgd
                   jc=pawfgrtab(iatom)%ifftsph(ic)
                   vloc(1,ic)=vtrial(jc,ispden)
                   vloc(2,ic)=zero
                 end do
               end if
             else if (ispden==3) then
               if (usexcnhat==0) then
                 vloc(:,:)=zero
               else
                 do ic=1,nfgd
                   jc=pawfgrtab(iatom)%ifftsph(ic)
                   vloc(1,ic)=vtrial(jc,3)
                   vloc(2,ic)=vtrial(jc,4)
                 end do
               end if
             else ! ispden=4
               vloc(2,1:nfgd)=-vloc(2,1:nfgd)
             end if
           end if

!          option = 0 Insulator case
           if(option==0)then
             do ilslm=1,lm_size
               do ic=1,nfgd
                 jc=pawfgrtab(iatom)%ifftsph(ic)
                 contrib(1:cplex_nspden) = zero

!                Int_R^3{vtrial*Sum_LM[Q_ij_q^LM^(1)]}
                 mua=alpha(istr);mub=beta(istr)
                 contrib(1:cplex_nspden)=contrib(1:cplex_nspden)+half*vloc(1:cplex_nspden,ic)&
&                  *(pawfgrtab(iatom)%gylmgr(mua,ic,ilslm)*pawfgrtab(iatom)%rfgd(mub,ic)&
&                  + pawfgrtab(iatom)%gylmgr(mub,ic,ilslm)*pawfgrtab(iatom)%rfgd(mua,ic))

!                Int_R^3{Vloc^(1)*Sum_LM[Q_ij_q^LM]}
                 contrib(1)=contrib(1)+vpsp1(jc)*pawfgrtab(iatom)%gylm(ic,ilslm)

!                delta_{alphabeta}Int_R^3{Vloc*Sum_LM[Q_ij_q^LM]}
                 if(istr<=3)then
                   contrib(1:cplex_nspden)=contrib(1:cplex_nspden) &
&                             +vloc(1:cplex_nspden,ic)*pawfgrtab(iatom)%gylm(ic,ilslm)
                 end if

                 intv(1:cplex_nspden,ilslm)=intv(1:cplex_nspden,ilslm)+contrib(1:cplex_nspden)
               end do
             end do

!          option = 1 Metal case (without Vpsp1)
           else if (option==1)then
             do ilslm=1,lm_size
               do ic=1,nfgd
                 jc=pawfgrtab(iatom)%ifftsph(ic)
                 contrib(1) = zero

!                Int_R^3{vtrial*Sum_LM[Q_ij_q^LM^(1)]}
                 mua=alpha(istr);mub=beta(istr)
                 contrib(1:cplex_nspden)=contrib(1:cplex_nspden)+half*vloc(1:cplex_nspden,ic)&
&                  *(pawfgrtab(iatom)%gylmgr(mua,ic,ilslm)*pawfgrtab(iatom)%rfgd(mub,ic)&
&                  + pawfgrtab(iatom)%gylmgr(mub,ic,ilslm)*pawfgrtab(iatom)%rfgd(mua,ic))

!                delta_{alphabeta}Int_R^3{Vtrial*Sum_LM[Q_ij_q^LM]}
                 if(istr<=3)then
                   contrib(1:cplex_nspden)=contrib(1:cplex_nspden) &
&                             +vloc(1:cplex_nspden,ic)*pawfgrtab(iatom)%gylm(ic,ilslm)
                 end if

                 intv(1:cplex_nspden,ilslm)=intv(1:cplex_nspden,ilslm)+contrib(1:cplex_nspden)
               end do
             end do
           end if
           LIBPAW_DEALLOCATE(vloc)

!          Apply ucvol/nfft factor on integral
           intv(:,:)=fact*intv(:,:)

!          --- Reduction in case of parallelization ---
           call xmpi_sum(intv,my_comm_grid,ier)

           paw_ij1(iatom)%dijfr(:,ispden)=zero

!          ---- Loop over (i,j) components
           klmn1=1
           do klmn=1,lmn2_size
             klm =pawtab(itypat)%indklmn(1,klmn)
             lmin=pawtab(itypat)%indklmn(3,klmn)
             lmax=pawtab(itypat)%indklmn(4,klmn)
             do ils=lmin,lmax,2
               lm0=ils**2+ils+1
               do mm=-ils,ils
                 ilslm=lm0+mm;isel=pawang%gntselect(lm0+mm,klm)
                 if (isel>0) then
                   !The following works only because cplex_nspden<=cplex_dij
                   paw_ij1(iatom)%dijfr(klmn1:klmn1+dplex_nsp,ispden)= &
&                    paw_ij1(iatom)%dijfr(klmn1:klmn1+dplex_nsp,ispden) &
&                    +pawtab(itypat)%qijl(ilslm,klmn)*intv(1:cplex_nspden,ilslm)
                 end if
               end do
             end do
             klmn1=klmn1+cplex_dij
           end do

!          Dijfr is marked as computed
           paw_ij1(iatom)%has_dijfr=2

         end if ! ipert

         LIBPAW_DEALLOCATE(intv)

!----------------------------------------------------------
!      End loops over spin components
       end do ! ispden

!      ----------------------------------------------------------
!      Deduce some part of Dij according to symmetries
!      ----------------------------------------------------------

     !Non-collinear: D_ij(:,4)=D^21_ij=D^12_ij^*
     else if (nspden==4.and.idij==4) then
       paw_ij1(iatom)%dijfr(:,idij)=paw_ij1(iatom)%dijfr(:,idij-1)
       if (cplex_dij==2) then
         do klmn=2,lmn2_size*cplex_dij,cplex_dij
           paw_ij1(iatom)%dijfr(klmn,idij)=-paw_ij1(iatom)%dijfr(klmn,idij)
         end do
         if (qphase==2) then
           do klmn=2+lmn2_size*cplex_dij,2*lmn2_size*cplex_dij,cplex_dij
             paw_ij1(iatom)%dijfr(klmn,idij)=-paw_ij1(iatom)%dijfr(klmn,idij)
           end do
         end if
       end if

     !Antiferro: D_ij(:,2)=D^down_ij=D^up_ij
     else if (nsppol==1.and.idij==2) then
       paw_ij1(iatom)%dijfr(:,idij)=paw_ij1(iatom)%dijfr(:,idij-1)
     end if

!  End loop on Dij components
   end do ! idij

!----------------------------------------------------------

!  Eventually free temporary space for g_l(r).Y_lm(r) gradients and exp(-i.q.r)
   if (need_dijfr_1.or.need_dijfr_2) then
     if (pawfgrtab(iatom)%gylm_allocated==2) then
       LIBPAW_DEALLOCATE(pawfgrtab(iatom)%gylm)
       LIBPAW_ALLOCATE(pawfgrtab(iatom)%gylm,(0,0))
       pawfgrtab(iatom)%gylm_allocated=0
     end if
     if (pawfgrtab(iatom)%gylmgr_allocated==2) then
       LIBPAW_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
       LIBPAW_ALLOCATE(pawfgrtab(iatom)%gylmgr,(0,0,0))
       pawfgrtab(iatom)%gylmgr_allocated=0
     end if
     if (pawfgrtab(iatom)%expiqr_allocated==2) then
       LIBPAW_DEALLOCATE(pawfgrtab(iatom)%expiqr)
       LIBPAW_ALLOCATE(pawfgrtab(iatom)%expiqr,(0,0))
       pawfgrtab(iatom)%expiqr_allocated=0
     end if
   end if

!  End loop on atoms
 end do

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine pawdijfr
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawpupot
!! NAME
!! pawpupot
!!
!! FUNCTION
!! Compute the PAW LDA+U on-site potential
!!
!! INPUTS
!!  cplex_dij=2 if LDA+U pot. is COMPLEX (as in the spin-orbit case), 1 if dij is REAL
!!  ndij=number of spin components for Dij
!!  pawprtvol=control print volume and debugging output for PAW
!!  noccmmp(cplex_dij,2*lpawu+1,2*lpawu+1,ndij)=density matrix in the augm. region
!!  nocctot(ndij)=number of electrons in the correlated subspace
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data:
!!     %usepawu, %upawu, %jpau
!!     %vee(2*lpawu+1*4)=screened coulomb matrix
!!
!! OUTPUT
!!  vpawu(cplex_dij,lpawu*2+1,lpawu*2+1,ndij)=lda+u potential
!!                                 (see eg PRB 52, 5467 (1995) [[cite:Liechenstein1995]])
!!    When vpawu is complex (cplex_dij=2):
!!      vpawu(2*i-1,:) contains the real part
!!      vpawu(2*i,:) contains the imaginary part
!!
!! PARENTS
!!      ldau_self,m_pawdij,m_pawhr
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

 subroutine pawpupot(cplex_dij,ndij,noccmmp,nocctot,&
&                    pawprtvol,pawtab,vpawu)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex_dij,ndij,pawprtvol
 type(pawtab_type),intent(in) :: pawtab
!arrays
 real(dp),intent(in) :: noccmmp(:,:,:,:),nocctot(:)
 real(dp),intent(out) :: vpawu(:,:,:,:)

!Local variables ---------------------------------------
!scalars
!Option for interaction energy in case of non-collinear magnetism:
!           1: E_int=-J/4.N.(N-2)                   (better)
!           2: E_int=-J/2.(Nup.(Nup-1)+Ndn.(Ndn-1)) (Nup and Ndn are ill-defined)
 integer,parameter :: option_interaction=1

 integer :: iplex,ispden,jspden,lpawu,m1,m11,m2,m21,m3,m31,m4,m41,nspden_eff
 real(dp) :: mnorm,mx,my,mz,n_sig,n_msig,n_tot,VUKStemp,n_sigs,n_msigs
 real(dp),save :: VUKS
 character(len=500) :: msg
!arrays
 real(dp),parameter :: factcg(3:4)=(/one,-one/)
 real(dp) :: n34_msig(cplex_dij),n34_sig(cplex_dij)
!real(dp) :: n43_sig(cplex_dij)

! *****************************************************

!Useful data
 lpawu=pawtab%lpawu

!Check data consistency
 if(pawtab%usepawu<0) then
   msg = "usepawu<0 not allowed!"
   MSG_BUG(msg)
 end if
 if(option_interaction==3.and.pawtab%usepawu>=10) then
   msg = "Option_interaction==3 is not compatible with usepawu>=10 in pawpupot"
   MSG_ERROR(msg)
 end if
 if (size(vpawu,1)/=cplex_dij.or.size(vpawu,2)/=2*lpawu+1.or.&
&    size(vpawu,3)/=2*lpawu+1.or.size(vpawu,4)/=ndij) then
   write (msg,'(a,4I5, a,4I5)') ' invalid sizes for vpawu !',cplex_dij,2*lpawu+1,2*lpawu+1,ndij, &
&    ' /= ', size(vpawu,1), size(vpawu,2), size(vpawu,3), size(vpawu,4)
   MSG_BUG(msg)
 end if
 if (size(noccmmp,1)/=cplex_dij.or.size(noccmmp,2)/=2*lpawu+1.or.&
&    size(noccmmp,3)/=2*lpawu+1.or.size(noccmmp,4)/=ndij) then
   write (msg,'(a,4I5, a,4I5)') ' invalid sizes for noccmmp !',cplex_dij,2*lpawu+1,2*lpawu+1,ndij, &
&    ' /= ', size(noccmmp,1), size(noccmmp,2), size(noccmmp,3), size(noccmmp,4)
   MSG_BUG(msg)
 end if
 if (size(nocctot,1)/=ndij) then
   msg='invalid size for nocctot !'
   MSG_BUG(msg)
 end if

!=====================================================
!Compute LDA+U Potential on the basis of projectors
!cf PRB 52 5467 (1995) [[cite:Liechenstein1995]]
!-----------------------------------------------------

 vpawu=zero ; nspden_eff=ndij
 do ispden=1,nspden_eff

   if (ispden<=2) then   ! cases ndij=4, ispden=1,2 or ndij<4
     jspden=min(nspden_eff,2)-ispden+1   ! (ispden,ndij)=(1,4)=>jspden=2

     if (nspden_eff<=2) then
       n_sig =nocctot(ispden)
       n_msig=nocctot(jspden)
       n_tot =n_sig+n_msig
     else
       n_tot=nocctot(1)
       mx=nocctot(2)
       my=nocctot(3)
       mz=nocctot(4)
       mnorm=sqrt(mx*mx+my*my+mz*mz)
       if (ispden==1) then
!        n_sig =half*(n_tot+mnorm)
!        n_msig=half*(n_tot-mnorm)
         n_sig =half*(n_tot+sign(mnorm,mz))
         n_msig=half*(n_tot-sign(mnorm,mz))
       else
!        n_sig =half*(n_tot-mnorm)
!        n_msig=half*(n_tot+mnorm)
         n_sig =half*(n_tot-sign(mnorm,mz))
         n_msig=half*(n_tot+sign(mnorm,mz))
       end if
     end if

     n_sigs =n_sig/(float(2*lpawu+1))
     n_msigs =n_msig/(float(2*lpawu+1))
     do m1=-lpawu,lpawu
       m11=m1+lpawu+1
       do m2=-lpawu,lpawu
         m21=m2+lpawu+1
         do m3=-lpawu,lpawu
           m31=m3+lpawu+1
           do m4=-lpawu,lpawu
             m41=m4+lpawu+1
             n34_sig(:) =noccmmp(:,m31,m41,ispden) ! spin sigma
             n34_msig(:)=noccmmp(:,m31,m41,jspden) ! opposite spin (-sigma)
             if(m31==m41.and.pawtab%usepawu==3) then
               n34_sig(1)= n34_sig(1) - n_sigs
               n34_msig(1)= n34_msig(1) - n_msigs
             end if
             do iplex=1,cplex_dij
               vpawu(iplex,m11,m21,ispden)=vpawu(iplex,m11,m21,ispden) &
&               +n34_msig(iplex)*pawtab%vee(m11,m31,m21,m41) &
&             +n34_sig(iplex)*(pawtab%vee(m11,m31,m21,m41)-pawtab%vee(m11,m31,m41,m21))
             end do
!            if(abs(pawprtvol)>=3.and.m11==1.and.m21==1) then
!            write(msg,'(a,i4,i4,2e20.10)') "m31,m41,vu=",m31,m41,&
!            & vpawu(:,m11,m21,ispden)
!            call wrtout(std_out,msg,'COLL')
!            write(msg,'(a,4e20.10)') "vee",pawtab%vee(m11,m31,m21,m41),&
!            & pawtab%vee(m11,m31,m41,m21)
!            call wrtout(std_out,msg,'COLL')
!            write(msg,'(a,4e20.10)') "n34_msig,n34_sig",n34_msig(1),n34_sig(1)
!            call wrtout(std_out,msg,'COLL')
!            end if
           end do
         end do
!        if(abs(pawprtvol)>=3) then
!        if(m11/=m21) then
!        write(msg,'(a,i4,i4,2e20.10)') "vu=",m11,m21,vpawu(:,m11,m21,ispden)
!        call wrtout(std_out,msg,'COLL')
!        write(msg,'(a,2e20.10)') "vupred=",-pawtab%upawu*noccmmp(:,m21,m11,ispden)
!        call wrtout(std_out,msg,'COLL')
!        end if
!        end if
       end do ! m2
       if(abs(pawprtvol)>=3) then
         write(msg,'(a,i3,14f11.5)') &
&         "vpawu   ",m11, (vpawu(:,m11,m21,ispden),m21=1,2*lpawu+1)
         call wrtout(std_out,  msg,'COLL')
         write(msg,'(a,i3,14f11.5)') &
&         "noccmmp ",m11, (noccmmp(:,m11,m21,ispden),m21=1,2*lpawu+1)
         call wrtout(std_out,  msg,'COLL')
       end if

!      Full localized limit
       if(pawtab%usepawu==1.or.pawtab%usepawu==4) then ! not activated if usepawu=10 !!
!        Here we compute vpawu=vpawu-v_dc
         vpawu(1,m11,m11,ispden)=vpawu(1,m11,m11,ispden)-pawtab%upawu*(n_tot-half)
         if (ndij/=4.or.option_interaction==2) then
           if(pawtab%usepawu/=4) then
             vpawu(1,m11,m11,ispden)=vpawu(1,m11,m11,ispden)+pawtab%jpawu*(n_sig-half)
           else
             vpawu(1,m11,m11,ispden)=vpawu(1,m11,m11,ispden)+half*pawtab%jpawu*(n_tot-one)
           endif
         else if (ndij==4.and.option_interaction==1) then
           vpawu(1,m11,m11,ispden)=vpawu(1,m11,m11,ispden)+half*pawtab%jpawu*(n_tot-one)
         else if (ndij==4.and.option_interaction==3) then
!          Here vdc^{alpha,beta}=\vect{m}.\vect{sigma}^{\beta,\alpha}
           vpawu(1,m11,m11,ispden)=vpawu(1,m11,m11,ispden)+half*pawtab%jpawu*(n_tot-one)
         end if

!        Around mean field
       else if(pawtab%usepawu==2) then
         vpawu(1,m11,m11,ispden)=vpawu(1,m11,m11,ispden)-n_msig*pawtab%upawu &
&         -n_sig*(pawtab%upawu-pawtab%jpawu) &
&         *(dble(2*lpawu)/dble(2*lpawu+1))
       end if

!      if (abs(pawprtvol)>=3) then
!      write(msg,'(a,i4,i4,2x,e20.10)') "vudiag= ",m11,m11,vpawu(1,m11,m11,ispden)
!      call wrtout(std_out,  msg,'COLL')
!      write(msg,'(a,2e20.10)') "vudiagpred= ",pawtab%upawu*(half-noccmmp(:,m11,m11,ispden))
!      call wrtout(std_out,  msg,'COLL')
!      end if
       if(abs(pawprtvol)>=3) then
         write(msg,*) "nocctot",nocctot
         call wrtout(std_out,  msg,'COLL')
         write(msg,'(a,i3,14f11.5)') &
&         "vpawu  2",m11, (vpawu(:,m11,m21,ispden),m21=1,2*lpawu+1)
         call wrtout(std_out,  msg,'COLL')
         write(msg,'(a,i3,14f11.5)') &
&         "noccmmp2",m11, (noccmmp(:,m11,m21,ispden),m21=1,2*lpawu+1)
         call wrtout(std_out,  msg,'COLL')
       end if
     end do ! m1

   end if ! ispden<=2

!  Non-collinear magnetism: add non-diagonal term; see (Eq 6) in PRB 72, 024458 (2005) [[cite:Shurikov2005]]
!  BA Here, we compute the transpose --- with respect to spin indices --- of
!  BA equation (6) of this reference, because of differences in notations,
!  BA namely Eband=\sum rhoij^{alpha,beta}*Dij(beta,alpha) contrary to PRB 72, 024458 (2005) [[cite:Shurikov2005]]
   if (ispden>=3) then
     mx=nocctot(2)
     my=nocctot(3)
     do m1=-lpawu,lpawu
       m11=m1+lpawu+1
       do m2=-lpawu,lpawu
         m21=m2+lpawu+1
         do m3=-lpawu,lpawu
           m31=m3+lpawu+1
           do m4=-lpawu,lpawu
             m41=m4+lpawu+1
!            n43_sig(:) =noccmmp(:,m41,m31,ispden)
!            vpawu(1,m11,m21,ispden)=vpawu(1,m11,m21,ispden)-n43_sig(1)*pawtab%vee(m11,m31,m41,m21)
!            vpawu(2,m11,m21,ispden)=vpawu(2,m11,m21,ispden)+n43_sig(2)*pawtab%vee(m11,m31,m41,m21)
             n34_sig(:) =noccmmp(:,m31,m41,ispden)
             vpawu(1,m11,m21,ispden)=vpawu(1,m11,m21,ispden)-n34_sig(1)*pawtab%vee(m11,m31,m41,m21)
             vpawu(2,m11,m21,ispden)=vpawu(2,m11,m21,ispden)-n34_sig(2)*pawtab%vee(m11,m31,m41,m21)
           end do
         end do
       end do
       if((pawtab%usepawu==1.or.pawtab%usepawu==4).and.option_interaction==3) then ! not activated if usepawu=10 !!
         vpawu(1,m11,m11,ispden)=vpawu(1,m11,m21,ispden)+half*pawtab%jpawu*mx
         if(ispden==3) then
           vpawu(2,m11,m11,ispden)=vpawu(1,m11,m21,ispden)-half*pawtab%jpawu*my
         else
           vpawu(2,m11,m11,ispden)=vpawu(1,m11,m21,ispden)+half*pawtab%jpawu*my
         end if
       end if
     end do
   end if

   if(abs(pawprtvol)>=3) then
     write(std_out,*) "vpawu, ispden",ispden
     do m11=1,2*lpawu+1
       write(msg,'(12(1x,9(1x,"(",f10.7,",",f10.7,")")))') &
&         (vpawu(1:cplex_dij,m11,m21,ispden),m21=1,2*lpawu+1)
       call wrtout(std_out,msg,'COLL')
     end do
   end if

!  Printing for test
   if (abs(pawprtvol)>=3) then
     if (ispden==1) VUKS=zero
     VUKStemp=zero
     do m1=-lpawu,lpawu
       m11=m1+lpawu+1
       do m2=-lpawu,lpawu
         m21=m2+lpawu+1
         VUKStemp=VUKStemp+vpawu(1,m11,m21,ispden)*noccmmp(1,m11,m21,ispden)
         if (cplex_dij == 2) then
           VUKStemp=VUKStemp-vpawu(2,m11,m21,ispden)*noccmmp(2,m11,m21,ispden)
         end if
         write(msg,'(a,2e20.10,2e20.10)') "m1,m2,vpawu,noccmmp= ", &
&          vpawu(:,m11,m21,ispden),noccmmp(:,m11,m21,ispden)
         call wrtout(std_out,  msg,'COLL')
       end do
     end do
     VUKS=VUKS+VUKStemp
     write(msg,*) "pawpupot: VUKStemp= ",ispden,VUKStemp
     call wrtout(std_out,  msg,'COLL')
     if (ispden==nspden_eff) then
       write(msg,*) "pawpupot: VUKS= ",ispden,VUKS
       call wrtout(std_out,  msg,'COLL')
     end if
   end if

 end do ! Loop on ispden

 end subroutine pawpupot
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawxpot
!! NAME
!! pawxpot
!!
!! FUNCTION
!! Compute the PAW Local Exact-Exchange on-site potential
!!
!! INPUTS
!!  ndij=number of spin components for Dij
!!  pawprtvol=control print volume and debugging output for PAW
!!  paw_ij <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawtab <type(pawtab_type)>=paw tabulated starting data:
!!  pawrhoij <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!
!! OUTPUT
!!  paw_ij%vpawx(pawtab%lexexch*2+1,pawtab%lexexch*2+1)=local exact-exchange potential
!!
!! PARENTS
!!      m_pawdij,pawdenpot
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

 subroutine pawxpot(ndij,pawprtvol,pawrhoij,pawtab,vpawx)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ndij,pawprtvol
 type(pawrhoij_type),intent(in) :: pawrhoij
 type(pawtab_type),intent(in) :: pawtab
 real(dp),intent(out) :: vpawx(:,:,:)

!Local variables ---------------------------------------
!scalars
 integer :: cplex_rhoij,irhoij,irhoij1,ispden,jrhoij,jrhoij1,klmn,klmn1,lexexch,ll,lmn2_size
 integer :: m11,m21,m31,m41,n1,n2,n3,n4,nk,nn1,nn2,nspden_eff
 real(dp) :: tot
 character(len=500) :: msg
!arrays
 integer :: indn(3,3)
 real(dp) :: factnk(6)

! *****************************************************

!Useful data
 lexexch=pawtab%lexexch
 cplex_rhoij=pawrhoij%cplex_rhoij
 lmn2_size=pawtab%lmn2_size
 if (pawtab%nproju==1) nk=1
 if (pawtab%nproju==2) nk=6
 factnk(1)=one;factnk(2)=one;factnk(3)=one
 factnk(4)=two;factnk(5)=two;factnk(6)=two
 indn(1,1)=1;indn(1,2)=4;indn(1,3)=5
 indn(2,1)=4;indn(2,2)=2;indn(2,3)=6
 indn(3,1)=5;indn(3,2)=6;indn(3,3)=3

!Check data consistency
 if (size(vpawx,1)/=1.or.size(vpawx,2)/=lmn2_size.or.&
&    size(vpawx,3)/=ndij) then
   msg='invalid sizes for vpawx !'
   MSG_BUG(msg)
 end if
 if (pawrhoij%qphase==2) then
   msg='pawxpot not compatible with qphase=2 (DFPT)!'
   MSG_BUG(msg)
 end if

!=====================================================
!Compute local exact exchange Potential
!on the basis of projectors.
!-----------------------------------------------------

 vpawx=zero ; nspden_eff=ndij
 do ispden=1,nspden_eff
   jrhoij=1
   do irhoij=1,pawrhoij%nrhoijsel
     klmn=pawrhoij%rhoijselect(irhoij)
     if(pawtab%indklmn(3,klmn)==0.and.pawtab%indklmn(4,klmn)==2*lexexch) then
       m11=pawtab%klmntomn(1,klmn);m21=pawtab%klmntomn(2,klmn)
       n1=pawtab%klmntomn(3,klmn);n2=pawtab%klmntomn(4,klmn)
       nn1=(n1*n2)/2+1
       jrhoij1=1
       do irhoij1=1,pawrhoij%nrhoijsel
         klmn1=pawrhoij%rhoijselect(irhoij1)
         if(pawtab%indklmn(3,klmn1)==0.and.pawtab%indklmn(4,klmn1)==2*lexexch) then
           m31=pawtab%klmntomn(1,klmn1);m41=pawtab%klmntomn(2,klmn1)
           n3=pawtab%klmntomn(3,klmn1);n4=pawtab%klmntomn(4,klmn1)
           nn2=(n3*n4)/2+1
           do ll=1,lexexch+1
             vpawx(1,klmn,ispden)=vpawx(1,klmn,ispden)&
&             -pawtab%vex(m11,m31,m41,m21,ll)*pawtab%dltij(klmn1) &
&             *pawtab%fk(indn(nn1,nn2),ll)*pawrhoij%rhoijp(jrhoij1,ispden)
           end do

         end if
         jrhoij1=jrhoij1+cplex_rhoij
       end do !irhoij1
     end if
     jrhoij=jrhoij+cplex_rhoij
   end do !irhoij
 end do !ispden

!Test
 if (abs(pawprtvol)>=2) then
   tot=zero
   do ispden=1,pawrhoij%nspden
     jrhoij=1
     do irhoij=1,pawrhoij%nrhoijsel
       klmn=pawrhoij%rhoijselect(irhoij)
       tot=tot+vpawx(1,klmn,ispden)*pawrhoij%rhoijp(jrhoij,ispden)*pawtab%dltij(klmn)
       jrhoij=jrhoij+cplex_rhoij
     end do
   end do
   write(msg, '(a,es22.15)' )" Vpawx: tot=",tot*half
   call wrtout(std_out,msg,'COLL')
 end if

 end subroutine pawxpot
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/symdij
!! NAME
!! symdij
!!
!! FUNCTION
!! Symmetrize PAW non-local strengths Dij
!! Symmetrize total Dij or one part of it
!!
!! INPUTS
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1).
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  ipert=index of perturbation if pawrhoij is a pertubed rhoij
!!        no meaning for ground-state calculations (should be 0)
!!  [mpi_atmtab(:)]=--optional-- indexes of the atoms treated by current proc
!!  [comm_atom]=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell
!!  nsym=number of symmetry elements in space group
!!  ntypat=number of types of atoms in unit cell.
!!  option_dij=choose which part of Dij has to be symmetrized (which paw_ij(:)%dijxxx):
!!             0: total dij (dij)
!!             1: dij due to compensation charge (dijhat)
!!             2: dij due to +U (dijU)
!!             3: dij XC (dijxc)
!!             4: dij XC due to compensation charge (dijxc_hat)
!!             5: dij XC valence only (dijxc_val)
!!             6: dij spin-orbit (dijso)
!!             7: dij exact exchange (dijexxc)
!!             8: dij, RF frozen part (dijfr)
!!             9: dij due to nuclear dipoles
!!             10: dij Hartree
!!             11: dij Fock
!!  paw_ij(natom)%cplex_dij=1 if dij are REAL, 2 if they are COMPLEX
!!  paw_ij(natom)%qphase=2 if exp^(-i.q.r) phase from RF at q<>0, 1 otherwise
!!  paw_ij(natom)%lmn_size=number of (l,m,n) elements for the paw basis
!!  paw_ij(natom)%nspden=number of spin-density components
!!  paw_ij(natom)%nsppol=number of independant spin-density components
!!  paw_ij(natom)%dij(lmn2_size,nspden)=non-symmetrized paw dij quantities
!!  pawang <type(pawang_type)>=angular mesh discretization and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  [qphon(3)]=--optional-- (RF calculations only) - wavevector of the phonon
!!  rprimd(3,3)=real space primitive translations.
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on
!!                   reciprocal space primitive translations
!!
!! SIDE EFFECTS
!!  paw_ij(natom)%dij???(cplex_dij*lmn2_size,nspden)=symmetrized dij quantities as output
!!
!! PARENTS
!!      bethe_salpeter,dfpt_scfcv,m_pawdij,paw_mknewh0,respfn,scfcv,screening
!!      sigma
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine symdij(gprimd,indsym,ipert,my_natom,natom,nsym,ntypat,option_dij,&
&                 paw_ij,pawang,pawprtvol,pawtab,rprimd,symafm,symrec, &
&                 mpi_atmtab,comm_atom,qphon) ! optional arguments (parallelism)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ipert,my_natom,natom,nsym,ntypat,option_dij,pawprtvol
 integer,optional,intent(in) :: comm_atom
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symafm(nsym),symrec(3,3,nsym)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3)
 real(dp),intent(in),optional :: qphon(3)
 type(paw_ij_type),intent(inout) :: paw_ij(my_natom)
 type(pawtab_type),target,intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: at_indx,cplex_dij,iafm,iatom,iatom_tot,ii
 integer :: il,il0,ilmn,iln,iln0,ilpm,indexi,indexii,indexj,indexjj,indexjj0,indexk,indexkc,indexkc_q
 integer :: iplex,iq,irot,ispden,itypat,j0lmn,jl,jl0,jlmn,jln,jln0,jlpm,jspden
 integer :: klmn,klmnc,kspden,lmn_size,lmn2_size,mi,mj,my_comm_atom,my_cplex_dij,my_ndij,my_qphase
 integer :: mu,natinc,ndij0,ndij1,nu,qphase,sz1,sz2
 logical,parameter :: afm_noncoll=.true.  ! TRUE if antiferro symmetries are used with non-collinear magnetism
 logical :: antiferro,has_qphase,my_atmtab_allocated,noncoll,paral_atom,use_afm
!DEBUG_ALTERNATE_ALGO
!Set to TRUE to choose an alternate algorithm (with another representation)
!to symmetrize Dij within non-collinear magnetism or spin-orbit
 logical,parameter :: lsymnew=.false.
!DEBUG_ALTERNATE_ALGO
 real(dp) :: arg,factafm,zarot2
 character(len=6) :: pertstrg,wrt_mode
 character(len=500) :: msg
!arrays
 integer :: nsym_used(2)
 integer, pointer :: indlmn(:,:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: dijc(2),fact(2),factsym(2),phase(2)
 real(dp) :: rotdij(2,2,2),rotmag(2,3,2),sumdij(2,2,2),summag(2,3,2)
 real(dp),allocatable :: dijnew(:,:,:),dijtmp(:,:),symrec_cart(:,:,:)
 type(coeff2_type),target, allocatable :: my_tmp_dij(:)
 type(coeff2_type),pointer :: tmp_dij(:)

!DEBUG_ALTERNATE_ALGO
!integer :: i1,i2,i3,i4,symrel_conv(3,3)
!real(dp) :: spinrot(4)
!real(dp),allocatable :: dijtemp(:,:),sumrhoso(:,:)
!complex(dpc) :: dijt(2,2),dijt2(2,2),Rspinrot(2,2)
!DEBUG_ALTERNATE_ALGO

! *********************************************************************

!Tests consistency of options
 if (my_natom>0) then
   if ((option_dij==1.and.paw_ij(1)%has_dijhat==0).or.&
&   (option_dij==2.and.paw_ij(1)%has_dijU==0).or.&
&   (option_dij==3.and.paw_ij(1)%has_dijxc==0).or.&
&   (option_dij==4.and.paw_ij(1)%has_dijxc_hat==0).or.&
&   (option_dij==5.and.paw_ij(1)%has_dijxc_val==0).or.&
&   (option_dij==6.and.paw_ij(1)%has_dijso==0).or.&
&   (option_dij==7.and.paw_ij(1)%has_dijexxc==0).or.&
&   (option_dij==8.and.paw_ij(1)%has_dijfr==0).or.&
&   (option_dij==9.and.paw_ij(1)%has_dijnd==0).or.&
&   (option_dij==10.and.paw_ij(1)%has_dijhartree==0).or.&
&   (option_dij==11.and.paw_ij(1)%has_dijfock==0)) then
     msg='Incompatibilty between option_dij and allocation of Dij!'
     MSG_BUG(msg)
   end if
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Determine sizes of dij array according to options
 my_qphase=1;my_cplex_dij=1;my_ndij=1
 if (my_natom>0) then
   my_qphase=paw_ij(1)%qphase
   my_cplex_dij=paw_ij(1)%cplex_dij
   my_ndij=paw_ij(1)%ndij
   if (option_dij==4.or.option_dij==5.or.option_dij==9) my_qphase=1
   if (option_dij==10) my_cplex_dij=1
   if (option_dij==10) my_ndij=1
 end if

!Antiferro case ?
 antiferro=.false.;if (my_natom>0) antiferro=(paw_ij(1)%nspden==2.and.paw_ij(1)%nsppol==1.and.my_ndij/=4)
! Non-collinear case
 noncoll=.false.;if (my_natom>0) noncoll=(my_ndij==4)
 if (my_natom>0) then
   if (noncoll.and.paw_ij(1)%cplex_dij/=2) then
     msg='cplex_dij must be 2 with ndij=4!'
     MSG_BUG(msg)
   end if
 end if
!Do we use antiferro symmetries ?
 use_afm=((antiferro).or.(noncoll.and.afm_noncoll))

!Do we have a phase due to q-vector?
 has_qphase=.false.
 if (my_natom>0) then
   has_qphase=(paw_ij(1)%qphase==2)
   if (present(qphon)) then
     if (any(abs(qphon(1:3))>tol8).and.(.not.has_qphase)) then
       msg='Should have qphase=2 for a non-zero q!'
       MSG_BUG(msg)
     end if
   end if
!DEBUG_ALTERNATE_ALGO
!  if(lsymnew.and.has_qphase) then
!    msg='symdij: alternate algo not available for phonons at q<>0!'
!    MSG_BUG(msg)
!  end if
!DEBUG_ALTERNATE_ALGO
 end if

!Printing of unsymetrized Dij
 if (abs(pawprtvol)>=1.and.option_dij==0.and.ipert/=natom+1.and.ipert/=natom+10) then
   wrt_mode='COLL';if (paral_atom) wrt_mode='PERS'
   pertstrg="DIJ";if (ipert>0) pertstrg="DIJ(1)"
   natinc=1;if(my_natom>1.and.pawprtvol>=0) natinc=my_natom-1
   write(msg, '(7a)') ch10," PAW TEST:",ch10,&
&     ' ========= Values of ',trim(pertstrg),' before symetrization =========',ch10
   call wrtout(std_out,msg,wrt_mode)
   do iatom=1,my_natom,natinc
     iatom_tot=iatom; if (paral_atom) iatom_tot=my_atmtab(iatom)
     call pawdij_print_dij(paw_ij(iatom)%dij,paw_ij(iatom)%cplex_dij,paw_ij(iatom)%qphase,&
&                iatom_tot,natom,paw_ij(iatom)%nspden,opt_prtvol=pawprtvol,mode_paral=wrt_mode)
   end do
   call wrtout(std_out,"",wrt_mode)
 end if

!Symmetrization occurs only when nsym>1
 if (nsym>1.and.ipert/=natom+1.and.ipert/=natom+10) then

   if (pawang%nsym==0) then
     msg='pawang%zarot must be allocated!'
     MSG_BUG(msg)
   end if

!  Have to make a temporary copy of dij
   LIBPAW_DATATYPE_ALLOCATE(my_tmp_dij,(my_natom))
   if (my_natom>0) then
     do iatom=1,my_natom
       lmn2_size=paw_ij(iatom)%lmn2_size
       sz1=my_qphase*my_cplex_dij*lmn2_size;sz2=my_ndij
       LIBPAW_ALLOCATE(my_tmp_dij(iatom)%value,(sz1,sz2))
       LIBPAW_ALLOCATE(dijtmp,(sz1,sz2))
       if (option_dij==0) then
         dijtmp(:,:)=paw_ij(iatom)%dij(:,:)
       else if (option_dij==1) then
         dijtmp(:,:)=paw_ij(iatom)%dijhat(:,:)
       else if (option_dij==2) then
         dijtmp(:,:)=paw_ij(iatom)%dijU(:,:)
       else if (option_dij==3) then
         dijtmp(:,:)=paw_ij(iatom)%dijxc(:,:)
       else if (option_dij==4) then
         dijtmp(:,:)=paw_ij(iatom)%dijxc_hat(:,:)
       else if (option_dij==5) then
         dijtmp(:,:)=paw_ij(iatom)%dijxc_val(:,:)
       else if (option_dij==6) then
         dijtmp(:,:)=paw_ij(iatom)%dijso(:,:)
       else if (option_dij==7) then
         dijtmp(:,:)=paw_ij(iatom)%dijexxc(:,:)
       else if (option_dij==8) then
         dijtmp(:,:)=paw_ij(iatom)%dijfr(:,:)
       else if (option_dij==9) then
         dijtmp(:,:)=paw_ij(iatom)%dijnd(:,:)
       else if (option_dij==10) then
         dijtmp(:,1)=paw_ij(iatom)%dijhartree(:)
       else if (option_dij==11) then
         dijtmp(:,:)=paw_ij(iatom)%dijfock(:,:)
       end if
       !Has to translate Dij^{alpha,beta} into (Dij, Dij magnetic field) format
       if (my_ndij==4) then
         my_tmp_dij(iatom)%value(:,1)=dijtmp(:,1)+dijtmp(:,2)
         my_tmp_dij(iatom)%value(:,2)=dijtmp(:,3)+dijtmp(:,4)
         my_tmp_dij(iatom)%value(:,4)=dijtmp(:,1)-dijtmp(:,2)
         do klmn=1,paw_ij(iatom)%lmn2_size
           my_tmp_dij(iatom)%value(2*klmn-1,3)=-dijtmp(2*klmn  ,3)+dijtmp(2*klmn  ,4)
           my_tmp_dij(iatom)%value(2*klmn  ,3)= dijtmp(2*klmn-1,3)-dijtmp(2*klmn-1,4)
         end do
!DEBUG_ALTERNATE_ALGO
!        if(lsymnew) my_tmp_dij(iatom)%value(:,:)=dijtmp(:,:)
!DEBUG_ALTERNATE_ALGO
       else
         my_tmp_dij(iatom)%value(:,:)=dijtmp(:,:)
       end if
       LIBPAW_DEALLOCATE(dijtmp)
     end do
   end if

!  Parallelism: gather all Dij
   if (paral_atom) then
     LIBPAW_DATATYPE_ALLOCATE(tmp_dij,(natom))
     call pawdij_gather(my_tmp_dij,tmp_dij,my_comm_atom,my_atmtab)
     do iatom=1,my_natom
       LIBPAW_DEALLOCATE(my_tmp_dij(iatom)%value)
     end do
     LIBPAW_DATATYPE_DEALLOCATE(my_tmp_dij)
   else
     tmp_dij=>my_tmp_dij
   end if

   if (noncoll) then
     LIBPAW_ALLOCATE(symrec_cart,(3,3,nsym))
     do irot=1,nsym
       symrec_cart(:,:,irot)=symdij_symcart(gprimd,rprimd,symrec(:,:,irot))
     end do
!DEBUG_ALTERNATE_ALGO
!    if(lsymnew) then
!      LIBPAW_ALLOCATE(sumrhoso,(my_cplex_dij,4))
!    end if
!DEBUG_ALTERNATE_ALGO
   end if

   ndij1=1
   if (antiferro) ndij1=2
   if (noncoll)   ndij1=4
   ndij1=min(ndij1,my_ndij)
   ndij0=ndij1-1
   LIBPAW_ALLOCATE(dijnew,(my_cplex_dij,ndij1,my_qphase))

!  Loops over atoms and spin components
   do iatom=1,my_natom
     iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
     itypat=paw_ij(iatom)%itypat
     lmn_size=paw_ij(iatom)%lmn_size
     lmn2_size=paw_ij(iatom)%lmn2_size
     cplex_dij=min(paw_ij(iatom)%cplex_dij,my_cplex_dij)
     qphase=min(paw_ij(iatom)%qphase,my_qphase)
     indlmn => pawtab(itypat)%indlmn

!DEBUG_ALTERNATE_ALGO
!    if (noncoll.and.lsymnew) then
!      LIBPAW_ALLOCATE(dijtemp,(cplex_dij,my_ndij))
!    end if
!DEBUG_ALTERNATE_ALGO

     do ispden=1,paw_ij(iatom)%nsppol
       jspden=min(3-ispden,paw_ij(iatom)%nsppol)

!      Loops over (il,im) and (jl,jm)
       jl0=-1;jln0=-1;indexj=1
       do jlmn=1,lmn_size
         jl=indlmn(1,jlmn)
         jlpm=1+jl+indlmn(2,jlmn)
         jln=indlmn(5,jlmn)
         if (jln/=jln0) indexj=indexj+2*jl0+1
         j0lmn=jlmn*(jlmn-1)/2
         il0=-1;iln0=-1;indexi=1
         do ilmn=1,jlmn
           il=indlmn(1,ilmn)
           ilpm=1+il+indlmn(2,ilmn)
           iln=indlmn(5,ilmn)
           if (iln/=iln0) indexi=indexi+2*il0+1
           klmn=j0lmn+ilmn;klmnc=cplex_dij*(klmn-1)

           nsym_used(:)=0

           rotdij(:,:,:)=zero
           if (noncoll) rotmag(:,:,:)=zero
!DEBUG_ALTERNATE_ALGO
!          if (noncoll.and.lsymnew) sumrhoso(:,:)=zero
!DEBUG_ALTERNATE_ALGO

!          Loop over symmetries
           do irot=1,nsym
!DEBUG_ALTERNATE_ALGO
!            if(lsymnew) then
!              call mati3inv(symrec(:,:,irot),symrel_conv)
!              call getspinrot(rprimd,spinrot,symrel_conv)
!              Rspinrot(1,1)=cmplx(spinrot(1),-spinrot(4))
!              Rspinrot(1,2)=cmplx(-spinrot(3),-spinrot(2))
!              Rspinrot(2,1)=cmplx(spinrot(3),-spinrot(2))
!              Rspinrot(2,2)=cmplx(spinrot(1),spinrot(4))
!            end if
!DEBUG_ALTERNATE_ALGO
             if ((symafm(irot)/=1).and.(.not.use_afm)) cycle
             kspden=ispden;if (symafm(irot)==-1) kspden=jspden
             iafm=1;if ((antiferro).and.(symafm(irot)==-1)) iafm=2
             factafm=dble(symafm(irot))

             nsym_used(iafm)=nsym_used(iafm)+1
             at_indx=indsym(4,irot,iatom_tot)

             if (has_qphase) then
               arg=two_pi*(qphon(1)*indsym(1,irot,iatom)+qphon(2)*indsym(2,irot,iatom) &
&                         +qphon(3)*indsym(3,irot,iatom))
               phase(1)=cos(arg);phase(2)=sin(arg)
             end if

             sumdij(:,:,:)=zero
             if (noncoll) summag(:,:,:)=zero

!            Accumulate values over (mi,mj) and symmetries
             do mj=1,2*jl+1
               indexjj=indexj+mj;indexjj0=indexjj*(indexjj-1)/2
               do mi=1,2*il+1
                 indexii=indexi+mi
                 factsym(:)=one
                 if (indexii<=indexjj) then
                   indexk=indexjj0+indexii
                   factsym(2)=one
                 else
                   indexk=indexii*(indexii-1)/2+indexjj
                   factsym(2)=-one
                 end if
                 indexkc=cplex_dij*(indexk-1)
                 indexkc_q=indexkc+cplex_dij*lmn2_size

!DEBUG_ALTERNATE_ALGO
!                if (noncoll.and.lsymnew) then
!                  do iplex=1,cplex_dij
!                    if(factafm>zero) then
!                      dijtemp(iplex,1)=factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,1)
!                      dijtemp(iplex,2)=factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,2)
!                    else
!                      dijtemp(iplex,1)=factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,2)
!                      dijtemp(iplex,2)=factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,1)
!                    end if
!                    if(factsym(2)<zero) then ! to be changed if symafm
!                      dijtemp(iplex,3)=factafm*factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,4)
!                      dijtemp(iplex,4)=factafm*factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,3)
!                    else
!                      dijtemp(iplex,3)=factafm*factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,3)
!                      dijtemp(iplex,4)=factafm*factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,4)
!                    end if
!                  end do
!                end if
!DEBUG_ALTERNATE_ALGO

!                Be careful: use here R_rel^-1 in term of spherical harmonics
!                which is tR_rec in term of spherical harmonics
!                so, use transpose[zarot]....  however, we use here zarot (??)
                 zarot2=pawang%zarot(mi,ilpm,il+1,irot)*pawang%zarot(mj,jlpm,jl+1,irot)
!                zarot2=pawang%zarot(ilpm,mi,il+1,irot)*pawang%zarot(jlpm,mj,jl+1,irot)

                 if((.not.noncoll).or.(.not.lsymnew)) then
                   fact(1)=factsym(1);fact(2)=factsym(2)*factafm   !????? What?  MT
                   sumdij(1:cplex_dij,iafm,1)=sumdij(1:cplex_dij,iafm,1) &
&                           +fact(1:cplex_dij)*zarot2 &
&                           *tmp_dij(at_indx)%value(indexkc+1:indexkc+cplex_dij,kspden)
                   if (qphase==2) &
&                    sumdij(1:cplex_dij,iafm,2)=sumdij(1:cplex_dij,iafm,2) &
&                             +fact(1:cplex_dij)*zarot2 &
&                             *tmp_dij(at_indx)%value(indexkc_q+1:indexkc_q+cplex_dij,kspden)
                 end if

                 if (noncoll.and.(.not.lsymnew)) then
                   fact(1)=factsym(1)*factafm;fact(2)=factsym(2)
                   do mu=1,3
                     summag(1:cplex_dij,mu,1)=summag(1:cplex_dij,mu,1) &
&                             +fact(1:cplex_dij)*zarot2 &
&                             *tmp_dij(at_indx)%value(indexkc+1:indexkc+cplex_dij,1+mu)
                   end do
                   if (qphase==2) then
                     do mu=1,3
                       summag(1:cplex_dij,mu,2)=summag(1:cplex_dij,mu,2) &
&                               +fact(1:cplex_dij)*zarot2 &
&                               *tmp_dij(at_indx)%value(indexkc_q+1:indexkc_q+cplex_dij,1+mu)
                     end do
                   end if
                 end if
!DEBUG_ALTERNATE_ALGO
!                if (noncoll.and.(lsymnew)) then
!                  dijt(1,1)=cmplx(dijtemp(1,1),dijtemp(2,1))
!                  dijt(2,2)=cmplx(dijtemp(1,2),dijtemp(2,2))
!                  dijt(1,2)=cmplx(dijtemp(1,3),dijtemp(2,3))
!                  dijt(2,1)=cmplx(dijtemp(1,4),dijtemp(2,4))
!                  dijt2(:,:)=czero
!                  do i1=1,2
!                    do i4=1,2
!                      do i2=1,2
!                        do i3=1,2
!                          dijt2(i1,i4)=dijt2(i1,i4)+Rspinrot(i1,i2)*dijt(i2,i3)*conjg(Rspinrot(i4,i3))
!                        end do
!                      end do
!                    end do
!                  end do
!                  do mu=1,4
!                    if(mu==1) then
!                      i1=1;i4=1
!                    else if(mu==2) then
!                      i1=2;i4=2
!                    else if(mu==3) then
!                      i1=1;i4=2
!                    else if(mu==4) then
!                      i1=2;i4=1
!                    end if
!                    sumrhoso(1,mu)=sumrhoso(1,mu)+zarot2*real(dijt2(i1,i4))
!                    sumrhoso(2,mu)=sumrhoso(2,mu)+zarot2*imag(dijt2(i1,i4))
!                  end do
!                end if
               end do ! mi
             end do ! mj
!DEBUG_ALTERNATE_ALGO

!            Apply phase for phonons
             if (has_qphase) then
               !Remember, Dij is stored as follows:
               ! Dij=  [Dij(2klmn-1)+i.Dij(2klmn)]
               !    +i.[Dij(lnm2_size+2klmn-1)+i.Dij(lmn2_size+2klmn)]
               if((.not.noncoll).or.(.not.lsymnew)) then
                 do iplex=1,cplex_dij
                   dijc(1)=sumdij(iplex,iafm,1)
                   dijc(2)=sumdij(iplex,iafm,2)
                   sumdij(iplex,iafm,1)=phase(1)*dijc(1)-phase(2)*dijc(2)
                   sumdij(iplex,iafm,2)=phase(1)*dijc(2)+phase(2)*dijc(1)
                 end do
               end if
               if (noncoll.and.(.not.lsymnew)) then
                 do iplex=1,cplex_dij
                   do mu=1,3
                     dijc(1)=summag(iplex,mu,1)
                     dijc(2)=summag(iplex,mu,2)
                     summag(iplex,mu,1)=phase(1)*dijc(1)-phase(2)*dijc(2)
                     summag(iplex,mu,2)=phase(1)*dijc(2)+phase(2)*dijc(1)
                   end do
                 end do
               end if
!DEBUG_ALTERNATE_ALGO
!              if (noncoll.and.(lsymnew) then
!                do mu=1,4
!                  sumrhoso(1,mu)=phase(1)*sumrhoso(1,mu)-phase(2)*sumrhoso(2,mu)
!                  sumrhoso(2,mu)=phase(1)*sumrhoso(2,mu)+phase(2)*sumrhoso(1,mu)
!                  end do
!                end do
!              end if
!DEBUG_ALTERNATE_ALGO
             end if

!            Add contribution of this rotation
             do iq=1,qphase
               rotdij(1:cplex_dij,iafm,iq)=rotdij(1:cplex_dij,iafm,iq) &
&                                         +sumdij(1:cplex_dij,iafm,iq)
             end do
             if (noncoll.and.(.not.lsymnew)) then
!              If non-collinear case, rotate Dij magnetization
!              Should use symrel^1 but use transpose[symrec] instead
               do iq=1,qphase
                 do nu=1,3
                   do mu=1,3
                     !We need the transpose ?
                     rotmag(1:cplex_dij,mu,iq)=rotmag(1:cplex_dij,mu,iq) &
&                       +symrec_cart(mu,nu,irot)*summag(1:cplex_dij,nu,iq)
                   end do
                 end do
               end do
             end if

           end do ! End loop over symmetries

           if((.not.noncoll).or.(.not.lsymnew)) then
!            Store new value of dij
             do iq=1,qphase
               do iplex=1,cplex_dij
                 dijnew(iplex,1,iq)=rotdij(iplex,1,iq)/nsym_used(1)
                 if (abs(dijnew(iplex,1,iq))<=tol10) dijnew(iplex,1,iq)=zero
               end do
             end do

!            Antiferromagnetic case: has to fill up "down" component of dij
             if (antiferro.and.nsym_used(2)>0) then
               do iq=1,qphase
                 do iplex=1,cplex_dij
                   dijnew(iplex,2,iq)=rotdij(iplex,2,iq)/nsym_used(2)
                   if (abs(dijnew(iplex,2,iq))<=tol10) dijnew(iplex,2,iq)=zero
                 end do
               end do
             end if
!DEBUG_ALTERNATE_ALGO
!          else if (noncoll.and.(lsymnew)) then
!            do mu=1,4
!              do iplex=1,cplex_dij
!                dijnew(iplex,mu,1)=sumrhoso(iplex,mu)/nsym_used(1)
!                if (abs(dijnew(iplex,mu,1))<=tol10) dijnew(iplex,mu,1)=zero
!              end do
!            end do
!DEBUG_ALTERNATE_ALGO
           end if

!          Non-collinear case: store new values of Dij magnetization
           if (noncoll.and.(.not.lsymnew)) then
!            Select on-zero elements
             do iq=1,qphase
               do mu=1,3
                 do iplex=1,cplex_dij
                   rotmag(iplex,mu,iq)=rotmag(iplex,mu,iq)/nsym_used(1)
                   if (abs(rotmag(iplex,mu,iq))<=tol10) rotmag(iplex,mu,iq)=zero
                 end do
               end do
             end do
!            Transfer back to Dij^{alpha,beta}
             if(.not.lsymnew) then
               !Remember: cplex_dij is 2 in that case
               do iq=1,qphase
                 dijnew(1,1,iq)=half*(dijnew(1,1,iq)+rotmag(1,3,iq))
                 dijnew(2,1,iq)=half*(dijnew(2,1,iq)+rotmag(2,3,iq))
                 dijnew(1,2,iq)=      dijnew(1,1,iq)-rotmag(1,3,iq)
                 dijnew(2,2,iq)=      dijnew(2,1,iq)-rotmag(2,3,iq)
                 dijnew(1,3,iq)=half*(rotmag(1,1,iq)+rotmag(2,2,iq))
                 dijnew(2,3,iq)=half*(rotmag(2,1,iq)-rotmag(1,2,iq))
                 dijnew(1,4,iq)=half*(rotmag(1,1,iq)-rotmag(2,2,iq))
                 dijnew(2,4,iq)=half*(rotmag(2,1,iq)+rotmag(1,2,iq))
               end do
             end if
           end if
!          Transfer new value of Dij in suitable pointer
           ii=klmnc
           do iq=1,qphase
             if (option_dij==0) then
               paw_ij(iatom)%dij(ii+1:ii+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1,iq)
             else if (option_dij==1) then
               paw_ij(iatom)%dijhat(ii+1:ii+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1,iq)
             else if (option_dij==2) then
               paw_ij(iatom)%dijU(ii+1:ii+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1,iq)
             else if (option_dij==3) then
               paw_ij(iatom)%dijxc(ii+1:ii+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1,iq)
             else if (option_dij==4) then
               paw_ij(iatom)%dijxc_hat(ii+1:ii+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1,iq)
             else if (option_dij==5) then
               paw_ij(iatom)%dijxc_val(ii+1:ii+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1,iq)
             else if (option_dij==6) then
               paw_ij(iatom)%dijso(ii+1:ii+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1,iq)
             else if (option_dij==7) then
               paw_ij(iatom)%dijexxc(ii+1:ii+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1,iq)
             else if (option_dij==8) then
               paw_ij(iatom)%dijfr(ii+1:ii+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1,iq)
             else if (option_dij==9) then
               paw_ij(iatom)%dijnd(ii+1:ii+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1,iq)
             else if (option_dij==10) then
               paw_ij(iatom)%dijhartree(ii+1:ii+cplex_dij)=dijnew(1:cplex_dij,1,iq)
             else if (option_dij==11) then
               paw_ij(iatom)%dijfock(ii+1:ii+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1,iq)
             end if
             ii=ii+lmn2_size*cplex_dij
           end do

           il0=il;iln0=iln  ! End loops over (il,im) and (jl,jm)
         end do
         jl0=jl;jln0=jln
       end do

     end do ! ispden

!DEBUG_ALTERNATE_ALGO
!    if (noncoll.and.lsymnew) then
!      LIBPAW_DEALLOCATE(dijtemp)
!    end if
!DEBUG_ALTERNATE_ALGO

   end do ! iatom

   LIBPAW_DEALLOCATE(dijnew)
   if (noncoll)  then
     LIBPAW_DEALLOCATE(symrec_cart)
!DEBUG_ALTERNATE_ALGO
!    if (lsymnew) then
!      LIBPAW_DEALLOCATE(sumrhoso)
!    end if
!DEBUG_ALTERNATE_ALGO
   end if

   if (paral_atom) then
     do iatom=1,natom
       LIBPAW_DEALLOCATE(tmp_dij(iatom)%value)
     end do
     LIBPAW_DATATYPE_DEALLOCATE(tmp_dij)
   else
     do iatom=1,my_natom
       LIBPAW_DEALLOCATE(my_tmp_dij(iatom)%value)
     end do
     LIBPAW_DATATYPE_DEALLOCATE(my_tmp_dij)
   end if

 else if (ipert/=natom+1.and.ipert/=natom+10) then  ! nsym>1

!  *********************************************************************
!  If nsym==1, only cut small components of dij

   if (antiferro) then
     msg='In the antiferromagnetic case, nsym cannot be 1'
     MSG_BUG(msg)
   end if

   do iatom=1,my_natom
     do ispden=1,my_ndij
       qphase=paw_ij(iatom)%qphase
       cplex_dij=paw_ij(iatom)%cplex_dij
       lmn2_size=paw_ij(iatom)%lmn2_size
       if (option_dij==0) then
         do klmn=1,lmn2_size*cplex_dij*qphase
           if (abs(paw_ij(iatom)%dij(klmn,ispden))<=tol10) paw_ij(iatom)%dij(klmn,ispden)=zero
         end do
       else if (option_dij==1) then
         do klmn=1,lmn2_size*cplex_dij*qphase
           if (abs(paw_ij(iatom)%dijhat(klmn,ispden))<=tol10) paw_ij(iatom)%dijhat(klmn,ispden)=zero
         end do
       else if (option_dij==2) then
         do klmn=1,lmn2_size*cplex_dij*qphase
           if (abs(paw_ij(iatom)%dijU(klmn,ispden))<=tol10) paw_ij(iatom)%dijU(klmn,ispden)=zero
         end do
       else if (option_dij==3) then
         do klmn=1,lmn2_size*cplex_dij*qphase
           if (abs(paw_ij(iatom)%dijxc(klmn,ispden))<=tol10) paw_ij(iatom)%dijxc(klmn,ispden)=zero
         end do
       else if (option_dij==4) then
         do klmn=1,lmn2_size*cplex_dij
           if (abs(paw_ij(iatom)%dijxc_hat(klmn,ispden))<=tol10) paw_ij(iatom)%dijxc_hat(klmn,ispden)=zero
         end do
       else if (option_dij==5) then
         do klmn=1,lmn2_size*cplex_dij
           if (abs(paw_ij(iatom)%dijxc_val(klmn,ispden))<=tol10) paw_ij(iatom)%dijxc_val(klmn,ispden)=zero
         end do
       else if (option_dij==6) then
         do klmn=1,lmn2_size*cplex_dij*qphase
           if (abs(paw_ij(iatom)%dijso(klmn,ispden))<=tol10) paw_ij(iatom)%dijso(klmn,ispden)=zero
         end do
       else if (option_dij==7) then
         do klmn=1,lmn2_size*cplex_dij*qphase
           if (abs(paw_ij(iatom)%dijexxc(klmn,ispden))<=tol10) paw_ij(iatom)%dijexxc(klmn,ispden)=zero
         end do
       else if (option_dij==8) then
         do klmn=1,lmn2_size*cplex_dij*qphase
           if (abs(paw_ij(iatom)%dijfr(klmn,ispden))<=tol10) paw_ij(iatom)%dijfr(klmn,ispden)=zero
         end do
       else if (option_dij==9) then
         do klmn=1,lmn2_size*cplex_dij
           if (abs(paw_ij(iatom)%dijnd(klmn,ispden))<=tol10) paw_ij(iatom)%dijnd(klmn,ispden)=zero
         end do
       else if (option_dij==10.and.ispden==1) then
         do klmn=1,lmn2_size*qphase
           if (abs(paw_ij(iatom)%dijhartree(klmn))<=tol10) paw_ij(iatom)%dijhartree(klmn)=zero
         end do
       else if (option_dij==11) then
         do klmn=1,lmn2_size*cplex_dij*qphase
           if (abs(paw_ij(iatom)%dijfock(klmn,ispden))<=tol10) paw_ij(iatom)%dijfock(klmn,ispden)=zero
         end do
       end if
     end do
   end do

 end if  ! nsym>1

!*********************************************************************
!Printing of Dij

 if (abs(pawprtvol)>=1.and.option_dij==0.and.ipert/=natom+1.and.ipert/=natom+10) then
   wrt_mode='COLL';if (paral_atom) wrt_mode='PERS'
   pertstrg="DIJ";if (ipert>0) pertstrg="DIJ(1)"
   natinc=1;if(my_natom>1.and.pawprtvol>=0) natinc=my_natom-1
   write(msg, '(7a)') ch10," PAW TEST:",ch10,&
&     ' ========= Values of ',trim(pertstrg),' after symetrization =========',ch10
   call wrtout(std_out,msg,wrt_mode)
   do iatom=1,my_natom,natinc
     iatom_tot=iatom; if (paral_atom) iatom_tot=my_atmtab(iatom)
     call pawdij_print_dij(paw_ij(iatom)%dij,paw_ij(iatom)%cplex_dij,paw_ij(iatom)%qphase,&
&                iatom_tot,natom,paw_ij(iatom)%nspden,opt_prtvol=pawprtvol,mode_paral=wrt_mode)
   end do
   call wrtout(std_out,"",wrt_mode)
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

!*********************************************************************
!Small function: convert a symmetry operation
!from reduced coordinates (integers) to cartesian coordinates (reals)
 contains
   function symdij_symcart(aprim,bprim,symred)

   real(dp) :: symdij_symcart(3,3)
   integer,intent(in) :: symred(3,3)
   real(dp),intent(in) :: aprim(3,3),bprim(3,3)
   integer :: ii,jj,kk
   real(dp) :: tmp(3,3)
   symdij_symcart=zero;tmp=zero
   do kk=1,3
     do jj=1,3
       do ii=1,3
         tmp(ii,jj)=tmp(ii,jj)+bprim(ii,kk)*dble(symred(jj,kk))
       end do
     end do
   end do
   do kk=1,3
     do jj=1,3
       do ii=1,3
         symdij_symcart(ii,jj)=symdij_symcart(ii,jj)+aprim(ii,kk)*tmp(jj,kk)
       end do
     end do
   end do
   end function symdij_symcart

end subroutine symdij
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/symdij_all
!! NAME
!! symdij_all
!!
!! FUNCTION
!! Symmetrize all contributions to PAW non-local strengths Dij
!!
!! INPUTS
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1).
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  ipert=index of perturbation if pawrhoij is a pertubed rhoij
!!        no meaning for ground-state calculations (should be 0)
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell
!!  nsym=number of symmetry elements in space group
!!  ntypat=number of types of atoms in unit cell.
!!  paw_ij(natom)%cplex_dij=1 if dij are REAL, 2 if they are COMPLEX
!!  paw_ij(natom)%qphase=2 if exp^(-i.q.r) phase from RF at q<>0, 1 otherwise
!!  paw_ij(natom)%lmn_size=number of (l,m,n) elements for the paw basis
!!  paw_ij(natom)%nspden=number of spin-density components
!!  paw_ij(natom)%nsppol=number of independant spin-density components
!!  paw_ij(natom)%dij(lmn2_size,nspden)=non-symmetrized paw dij quantities
!!  pawang <type(pawang_type)>=angular mesh discretization and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  rprimd(3,3)=real space primitive translations.
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on
!!                   reciprocal space primitive translations
!!
!! SIDE EFFECTS
!!  paw_ij(natom)%dij???(cplex_dij*qphase*lmn2_size,nspden)=symmetrized dij quantities as output
!!
!! PARENTS
!!      paw_mknewh0,screening,sigma
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine symdij_all(gprimd,indsym,ipert,my_natom,natom,nsym,ntypat,&
&                     paw_ij,pawang,pawprtvol,pawtab,rprimd,symafm,symrec,&
&                     mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ipert,my_natom,natom,nsym,ntypat,pawprtvol
 integer,optional,intent(in) :: comm_atom
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symafm(nsym),symrec(3,3,nsym)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3)
 type(paw_ij_type),intent(inout) :: paw_ij(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer,parameter :: MAX_NOPTS=12
 integer :: ii,option_dij,my_comm_atom,nopt
 logical :: my_atmtab_allocated,paral_atom
 character(len=500) :: msg
!arrays
 integer :: options(MAX_NOPTS)
 integer,pointer :: my_atmtab(:)

! *********************************************************************

 nopt = 0
 if (ANY(paw_ij(:)%has_dij==2)) then
   nopt = nopt + 1
   options(nopt) = 0
 end if

 if (ANY(paw_ij(:)%has_dijhat==2)) then
   nopt = nopt + 1
   options(nopt) = 1
 end if

 if (ANY(paw_ij(:)%has_dijU==2))   then
   nopt = nopt + 1
   options(nopt) = 2
 end if

 if (ANY(paw_ij(:)%has_dijxc==2)) then
   nopt = nopt + 1
   options(nopt) = 3
 end if

 if (ANY(paw_ij(:)%has_dijxc_hat==2)) then
   nopt = nopt + 1
   options(nopt) = 4
 end if

 if (ANY(paw_ij(:)%has_dijxc_val==2)) then
   nopt = nopt + 1
   options(nopt) = 5
 end if

 if (ANY(paw_ij(:)%has_dijso==2)) then
   nopt = nopt + 1
   options(nopt) = 6
 end if

 if (ANY(paw_ij(:)%has_dijexxc==2)) then
   nopt = nopt + 1
   options(nopt) = 7
 end if

 if (ANY(paw_ij(:)%has_dijfr==2)) then
   nopt = nopt + 1
   options(nopt) = 8
 end if

 if (ANY(paw_ij(:)%has_dijnd==2)) then
   nopt = nopt + 1
   options(nopt) = 9
 end if

 if (ANY(paw_ij(:)%has_dijhartree==2)) then
   nopt = nopt + 1
   options(nopt) = 10
 end if

 if (ANY(paw_ij(:)%has_dijfock==2)) then
   nopt = nopt + 1
   options(nopt) = 11
 end if

 if (ANY(paw_ij(:)%has_exexch_pot==2)) then
   nopt = nopt + 1
   options(nopt) = 10
   msg='symetrization of dij_exexch not coded!'
   MSG_ERROR(msg)
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

 do ii=1,nopt
   option_dij = options(ii)
   if (paral_atom) then
     call symdij(gprimd,indsym,ipert,my_natom,natom,nsym,ntypat,option_dij,&
&     paw_ij,pawang,pawprtvol,pawtab,rprimd,symafm,symrec,&
&     comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
   else
     call symdij(gprimd,indsym,ipert,my_natom,natom,nsym,ntypat,option_dij,&
&     paw_ij,pawang,pawprtvol,pawtab,rprimd,symafm,symrec)
   end if
 end do

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine symdij_all
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawdij_gather
!! NAME
!!  pawdij_gather
!!
!! FUNCTION
!!  Performs a ALLGATHER operation (over atomic sites) on Dij data
!!  stored as a 1D array of Dij arrays.
!!
!! INPUTS
!!  dij_in = coeff2d_type array containing the input Dij
!!  comm_atom= MPI communicator over atoms
!!  mpi_atmtab(:)=indexes of the atoms treated by current proc
!!
!! OUTPUT
!!  dij_out = coeff2d_type array containing the gathered Dij
!!
!! PARENTS
!!      m_pawdij
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdij_gather(dij_in,dij_out,comm_atom,mpi_atmtab)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm_atom
!arrays
 integer,intent(in) :: mpi_atmtab(:)
 type(coeff2_type),intent(in) :: dij_in(:)
 type(coeff2_type),intent(out) :: dij_out(:)

!Local variables-------------------------------
!scalars
 integer :: buf_dp_size,buf_dp_size_all,buf_int_size,buf_int_size_all
 integer :: dij_size,dij_size_out,ierr,ii,i2,indx_dp,indx_int,ival,n1,n2,nproc
!arrays
 integer :: bufsz(2)
 integer, allocatable :: buf_int(:),buf_int_all(:)
 integer, allocatable :: count_dp(:),count_int(:),count_tot(:),displ_dp(:),displ_int(:)
 integer, allocatable :: dimdij(:,:)
 real(dp),allocatable :: buf_dp(:),buf_dp_all(:)

! *************************************************************************

 nproc=xmpi_comm_size(comm_atom)
 dij_size=size(dij_in,dim=1)

 buf_dp_size=0
 LIBPAW_ALLOCATE(dimdij,(dij_size,2))
 do ii=1,dij_size
   dimdij(ii,1)=size(dij_in(ii)%value,dim=1)
   dimdij(ii,2)=size(dij_in(ii)%value,dim=2)
   buf_dp_size=buf_dp_size+dimdij(ii,1)*dimdij(ii,2)
 end do

!If only one proc, perform a single copy
 if (nproc==1) then
   do ii=1,dij_size
     ival=mpi_atmtab(ii)
     if (allocated(dij_out(ival)%value)) then
       LIBPAW_DEALLOCATE(dij_out(ival)%value)
     end if
     LIBPAW_ALLOCATE(dij_out(ival)%value,(n1,n2))
     dij_out(ii)%value=dij_in(ival)%value
   end do
   LIBPAW_DEALLOCATE(dimdij)
   return
 end if

!Fill in integer buffer
 buf_int_size=3*dij_size
 LIBPAW_ALLOCATE(buf_int,(buf_int_size))
 indx_int=1
 do ii=1,dij_size
   buf_int(indx_int  )=dimdij(ii,1)
   buf_int(indx_int+1)=dimdij(ii,2)
   buf_int(indx_int+2)=mpi_atmtab(ii)
   indx_int=indx_int+3
 end do

!Fill in real buffer
 LIBPAW_ALLOCATE(buf_dp,(buf_dp_size))
 indx_dp=1
 do ii=1,dij_size
   n1=dimdij(ii,1); n2=dimdij(ii,2)
   do i2=1,n2
     buf_dp(indx_dp:indx_dp+n1-1)=dij_in(ii)%value(1:n1,i2)
     indx_dp=indx_dp+n1
   end do
 end do

!Communicate (1 gather for integers, 1 gather for reals)
 LIBPAW_ALLOCATE(count_int,(nproc))
 LIBPAW_ALLOCATE(displ_int,(nproc))
 LIBPAW_ALLOCATE(count_dp ,(nproc))
 LIBPAW_ALLOCATE(displ_dp ,(nproc))
 LIBPAW_ALLOCATE(count_tot,(2*nproc))
 bufsz(1)=buf_int_size; bufsz(2)=buf_dp_size
 call xmpi_allgather(bufsz,2,count_tot,comm_atom,ierr)
 do ii=1,nproc
   count_int(ii)=count_tot(2*ii-1)
   count_dp (ii)=count_tot(2*ii)
 end do
 displ_int(1)=0;displ_dp(1)=0
 do ii=2,nproc
   displ_int(ii)=displ_int(ii-1)+count_int(ii-1)
   displ_dp (ii)=displ_dp (ii-1)+count_dp (ii-1)
 end do
 buf_int_size_all=sum(count_int)
 buf_dp_size_all =sum(count_dp)
 LIBPAW_DEALLOCATE(count_tot)
 LIBPAW_ALLOCATE(buf_int_all,(buf_int_size_all))
 LIBPAW_ALLOCATE(buf_dp_all ,(buf_dp_size_all))
 call xmpi_allgatherv(buf_int,buf_int_size,buf_int_all,count_int,displ_int,comm_atom,ierr)
 call xmpi_allgatherv(buf_dp ,buf_dp_size ,buf_dp_all ,count_dp ,displ_dp ,comm_atom,ierr)
 LIBPAW_DEALLOCATE(count_int)
 LIBPAW_DEALLOCATE(displ_int)
 LIBPAW_DEALLOCATE(count_dp)
 LIBPAW_DEALLOCATE(displ_dp)

!Retrieve gathered data
 dij_size_out=buf_int_size_all/3
 indx_int=1;indx_dp=1
 do ii=1,dij_size_out
   n1=buf_int_all(indx_int)
   n2=buf_int_all(indx_int+1)
   ival=buf_int_all(indx_int+2)
   indx_int=indx_int+3
   if (allocated(dij_out(ival)%value)) then
     LIBPAW_DEALLOCATE(dij_out(ival)%value)
   end if
   LIBPAW_ALLOCATE(dij_out(ival)%value,(n1,n2))
   do i2=1,n2
     dij_out(ival)%value(1:n1,i2)=buf_dp_all(indx_dp:indx_dp+n1-1)
     indx_dp=indx_dp+n1
   end do
 end do

 LIBPAW_DEALLOCATE(buf_dp_all)
 LIBPAW_DEALLOCATE(buf_int_all)
 LIBPAW_DEALLOCATE(buf_int)
 LIBPAW_DEALLOCATE(buf_dp)
 LIBPAW_DEALLOCATE(dimdij)

end subroutine pawdij_gather
!!***

!----------------------------------------------------------------------

!!****f* m_pawdij/pawdij_print_ij
!! NAME
!! pawdij_print_dij
!!
!! FUNCTION
!!  Print out the content of a Dij matrix (total Dij) in a suitable format
!!
!! INPUTS
!!  dij(cplex_dij*qphase*lmn2_size,ndij)= input matrix to be printed
!!  cplex_dij=1 if Dij is real, 2 if Dij is complex
!!  qphase=1 if Dij contains no RF phase, 2 if it contains a exp(-iqr) RF phase
!!  iatom=current atom
!!  natom=total number of atoms in the system
!!  nspden=number of spin density components
!!  [Ha_or_eV]= 1: output in hartrees, 2: output in eV
!!  [opt_prtvol]= >=0 if up to 12 components of _ij matrix have to be printed
!!                 <0 if all components of ij_ matrix have to be printed (optional)
!!  [mode_paral]= parallel printing mode (optional, default='COLL')
!!  [test_value]=(real number) if positive, print a warning when the magnitude of Dij is greater (optional)
!!  [title_msg]=message to print as title (optional)
!!  [unit]=the unit number for output (optional)
!!
!! OUTPUT
!! (Only writing)
!!
!! PARENTS
!!      m_pawdij
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawdij_print_dij(dij,cplex_dij,qphase,iatom,natom,nspden,&
&           test_value,title_msg,unit,Ha_or_eV,opt_prtvol,mode_paral) ! Optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex_dij,iatom,natom,nspden,qphase
 integer,optional,intent(in) :: Ha_or_eV,opt_prtvol,unit
 real(dp),intent(in),optional :: test_value
 character(len=4),optional,intent(in) :: mode_paral
 character(len=100),optional,intent(in) :: title_msg
!arrays
 real(dp),intent(in),target :: dij(:,:)

!Local variables-------------------------------
 character(len=7),parameter :: dspin(6)=(/"up     ","down   ","up-up  ","dwn-dwn","up-dwn ","dwn-up "/)
 integer :: idij,idij_sym,kk,lmn_size,lmn2_size,my_idij,my_idij_sym
 integer :: my_prtvol,my_unt,my_Ha_or_eV,ndij,tmp_cplex_dij
 real(dp) :: my_test_value,test_value_eff
 character(len=4) :: my_mode
 character(len=2000) :: msg
!arrays
 integer :: idum(0)
 real(dp),allocatable,target :: dij1(:),dij2(:)
 real(dp),pointer :: dij2p(:),dij2p_(:)

! *************************************************************************

!Optional arguments
 my_unt   =std_out ; if (PRESENT(unit      )) my_unt   =unit
 my_mode  ='COLL'  ; if (PRESENT(mode_paral)) my_mode  =mode_paral
 my_prtvol=1       ; if (PRESENT(opt_prtvol)) my_prtvol=opt_prtvol
 my_test_value=-one; if (PRESENT(test_value)) my_test_value=test_value
 my_Ha_or_eV=1     ; if (PRESENT(Ha_or_eV))   my_Ha_or_eV=Ha_or_eV

!Title
 if (present(title_msg)) then
   if (trim(title_msg)/='') then
     write(msg, '(2a)') ch10,trim(title_msg)
     call wrtout(my_unt,msg,my_mode)
   end if
 end if

!Inits
 ndij=size(dij,2)
 lmn2_size=size(dij,1)/(cplex_dij*qphase)
 lmn_size=int(dsqrt(two*dble(lmn2_size)))
 if (qphase==2) then
   LIBPAW_ALLOCATE(dij1,(2*lmn2_size))
   LIBPAW_ALLOCATE(dij2,(2*lmn2_size))
 end if

! === Loop over Dij components ===
 do idij=1,ndij

   idij_sym=idij;if (ndij==4.and.idij>2) idij_sym=7-idij

   !Subtitle
   if (natom>1.or.nspden>1.or.ndij==4) then
     if (nspden==1.and.ndij/=4) write(msg,'(a,i3)') ' Atom #',iatom
     if (nspden==2) write(msg,'(a,i3,a,i1)')' Atom #',iatom,' - Spin component ',idij
     if (ndij==4) write(msg,'(a,i3,2a)') ' Atom #',iatom,' - Component ',trim(dspin(idij+2*(ndij/4)))
     call wrtout(my_unt,msg,my_mode)
   end if

   !Select upper and lower triangular parts
   my_idij=min(size(dij,2),idij)
   my_idij_sym=min(size(dij,2),idij_sym)
   if (qphase==1) then
     tmp_cplex_dij=cplex_dij
     dij2p  => dij(1:cplex_dij*lmn2_size:1,my_idij)
     dij2p_ => dij(1:cplex_dij*lmn2_size:1,my_idij_sym)
   else
     tmp_cplex_dij=2
     if (cplex_dij==1) then
       do kk=1,lmn2_size
         dij1(2*kk-1)= dij(kk,my_idij)
         dij1(2*kk  )= dij(kk+lmn2_size,my_idij)
         dij2(2*kk-1)= dij(kk,my_idij_sym)
         dij2(2*kk  )=-dij(kk+lmn2_size,my_idij_sym)
       end do
     else
       do kk=1,lmn2_size
         dij1(2*kk-1)= dij(2*kk-1,my_idij)-dij(2*kk  +2*lmn2_size,my_idij)
         dij1(2*kk  )= dij(2*kk  ,my_idij)+dij(2*kk-1+2*lmn2_size,my_idij)
         dij2(2*kk-1)= dij(2*kk-1,my_idij_sym)+dij(2*kk  +2*lmn2_size,my_idij_sym)
         dij2(2*kk  )= dij(2*kk  ,my_idij_sym)-dij(2*kk-1+2*lmn2_size,my_idij_sym)
       end do
     end if
     dij2p => dij1 ; dij2p_ => dij2
   end if

   !Printing
    test_value_eff=-one;if(my_test_value>zero.and.idij==1) test_value_eff=my_test_value
    call pawio_print_ij(my_unt,dij2p,lmn2_size,tmp_cplex_dij,lmn_size,-1,idum,0,&
&                       my_prtvol,idum,test_value_eff,my_Ha_or_eV,&
&                       opt_sym=2,asym_ij=dij2p_,mode_paral=my_mode)

  end do !idij

 if (qphase==2) then
   LIBPAW_DEALLOCATE(dij1)
   LIBPAW_DEALLOCATE(dij2)
 end if

end subroutine pawdij_print_dij
!!***

!----------------------------------------------------------------------

END MODULE m_pawdij
!!***

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
!! Copyright (C) 2013-2017 ABINIT group (MT, FJ, BA, JWZ)
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
 public :: pawdij         ! Dij total
 public :: pawdijfock     ! Dij Fock exact-exchange
 public :: pawdijhartree  ! Dij Hartree
 public :: pawdijxc       ! Dij eXchange-Correlation (using (r,theta,phi) grid)
 public :: pawdijxcm      ! Dij eXchange-Correlation (using (l,m) moments)
 public :: pawdijhat      ! Dij^hat (compensation charge contribution)
 public :: pawdijnd       ! Dij nuclear dipole
 public :: pawdijso       ! Dij spin-orbit
 public :: pawdiju        ! Dij LDA+U
 public :: pawdijexxc     ! Dij local exact-exchange
 public :: pawdijfr       ! 1st-order frozen Dij
 public :: pawpupot       ! On-site LDA+U potential
 public :: pawxpot        ! On-site local exact-exchange potential
 public :: symdij         ! Symmetrize total Dij or one part of it
 public :: symdij_all     ! Symmetrize all contributions to Dij
 public :: pawdij_gather  ! Perform a allgather operation on Dij
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
!!  cplex=(RF calculations only) - 1 if RF 1st-order quantities are REAL, 2 if COMPLEX
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
!!  paw_ij(iatom)%dij(cplex_dij*lmn2_size,ndij)= total Dij terms (GS calculation, ipert=0)
!!                                               total 1st-order Dij terms (RF ccalc., ipert>0)
!!  May be complex if cplex_dij=2
!!        dij(:,:,1) contains Dij^up-up
!!        dij(:,:,2) contains Dij^dn-dn
!!        dij(:,:,3) contains Dij^up-dn (only if nspinor=2)
!!        dij(:,:,4) contains Dij^dn-up (only if nspinor=2)
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
&          mpi_atmtab,comm_atom,mpi_comm_grid,hyb_mixing,hyb_mixing_sr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdij'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,enunit,ipert,my_natom,natom,nfft,nfftot
 integer,intent(in) :: nspden,ntypat,pawprtvol,pawspnorb,pawxcdev
 integer,optional,intent(in) :: electronpositron_calctype
 integer,optional,intent(in) :: comm_atom,mpi_comm_grid,natvshift
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
 integer :: cplex_dij,iatom,iatom_tot,idij,ipositron,itypat,klmn,klmn1,lm_size,lmn2_size
 integer :: lpawu,my_comm_atom,my_comm_grid,natvshift_,ndij,nsploop,nsppol,usexcnhat
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
   if (paw_ij(1)%cplex/=paw_an(1)%cplex) then
     msg='paw_ij()%cplex and paw_an()%cplex must be equal!'
     MSG_BUG(msg)
   end if
   if (ipert<=0.and.paw_ij(1)%cplex/=1) then
     msg='cplex must be 1 for GS calculations!'
     MSG_BUG(msg)
   end if
   if (paw_ij(1)%cplex_dij<cplex) then
     msg='cplex_dij must be >= cplex!'
     MSG_BUG(msg)
   end if
   if (paw_ij(1)%cplex/=cplex) then
     msg='paw_ij()%cplex must be equal to cplex!'
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
   lm_size=paw_an(iatom)%lm_size
   lmn2_size=paw_ij(iatom)%lmn2_size
   ndij=paw_ij(iatom)%ndij
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
   dijnd_available=.false. ; dijnd_prereq=.true.
   if (has_nucdipmom) dijnd_available=(ipert<=0.and.any(abs(nucdipmom(:,iatom))>tol8))
!  DijSO: not available for RF, positron; only for spin-orbit ; VHartree and Vxc needed
   dijso_available=(pawspnorb>0.and.ipert<=0.and.ipositron/=1)
   dijso_prereq=(paw_ij(iatom)%has_dijso==2.or.&
&               (paw_an(iatom)%has_vhartree>0.and.paw_an(iatom)%has_vxc>0))
!  DijU: not available for RF, positron; only for LDA+U
   dijU_available=(pawtab(itypat)%usepawu>0.and.ipert<=0.and.ipositron/=1)
   dijU_prereq=(paw_ij(iatom)%has_dijU==2.or.paw_ij(iatom)%has_pawu_occ>0)
!  DijExxc: not available for RF, positron; only for local exact exch. ; Vxc_ex needed
   dijexxc_available=(pawtab(itypat)%useexexch>0.and.ipert<=0.and.ipositron/=1)
   dijexxc_prereq=(paw_ij(iatom)%has_dijexxc==2.or.paw_ij(iatom)%has_exexch_pot>0)
!  DijXC^hat: no condition ; Vxc needed
   dijxchat_available=.true.
   dijxchat_prereq=(paw_ij(iatom)%has_dijxc_hat==2.or.paw_an(iatom)%has_vxc>0)
!  DijXC_val: no condition ; Vxc_val needed
   dijxcval_available=.true.
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
       LIBPAW_ALLOCATE(paw_ij(iatom)%dij,(cplex_dij*lmn2_size,ndij))
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
       LIBPAW_ALLOCATE(paw_ij(iatom)%dijhartree,(cplex*lmn2_size))
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
       LIBPAW_ALLOCATE(paw_ij(iatom)%dijxc,(cplex_dij*lmn2_size,ndij))
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
       LIBPAW_ALLOCATE(paw_ij(iatom)%dijhat,(cplex_dij*lmn2_size,ndij))
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
       LIBPAW_ALLOCATE(paw_ij(iatom)%dijso,(cplex_dij*lmn2_size,ndij))
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
       LIBPAW_ALLOCATE(paw_ij(iatom)%dijU,(cplex_dij*lmn2_size,ndij))
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
!      LIBPAW_ALLOCATE(paw_ij(iatom)%dijxc_hat,(cplex_dij*lmn2_size,ndij))
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
!      LIBPAW_ALLOCATE(paw_ij(iatom)%dijxc_val,(cplex_dij*lmn2_size,ndij))
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
       if (dij_need) paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)+paw_ij(iatom)%dijfock(:,:)

     else

!    ===== Need to compute DijFock
       LIBPAW_ALLOCATE(dijfock_vv,(cplex_dij*lmn2_size,ndij))
       LIBPAW_ALLOCATE(dijfock_cv,(cplex_dij*lmn2_size,ndij))
       dijfock_vv(:,:)=zero ; dijfock_cv(:,:)=zero
!      Exact exchange is evaluated for electrons only
       if (ipositron/=1) then
         call pawdijfock(cplex,cplex_dij,dijfock_vv,dijfock_cv,hyb_mixing_,hyb_mixing_sr_, &
&                        ndij,nspden,nsppol,pawrhoij(iatom),pawtab(itypat))
       end if
       if (dijfock_need) paw_ij(iatom)%dijfock(:,:)=dijfock_vv(:,:)+dijfock_cv(:,:)
       if (dij_need) paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)+dijfock_vv(:,:)+dijfock_cv(:,:)
       LIBPAW_DEALLOCATE(dijfock_vv)
       LIBPAW_DEALLOCATE(dijfock_cv)
     end if
   end if

!  ----------- Add Dij_Hartree to Dij
!  ------------------------------------------------------------------------

   if ((dijhartree_need.or.dij_need).and.dijhartree_available) then

     LIBPAW_ALLOCATE(dijhartree,(cplex*lmn2_size))
!    ===== DijHartree already computed
     if (paw_ij(iatom)%has_dijhartree==2) then
       dijhartree(:)=paw_ij(iatom)%dijhartree(:)
     else
!    ===== Need to compute DijHartree
       if (ipositron/=1) then
         call pawdijhartree(cplex,dijhartree,nspden,pawrhoij(iatom),pawtab(itypat))
       else
         dijhartree(:)=zero
       end if
       if (ipositron/=0) then
         LIBPAW_ALLOCATE(dij_ep,(cplex*lmn2_size))
         call pawdijhartree(cplex,dij_ep,nspden,electronpositron_pawrhoij(iatom),pawtab(itypat))
         dijhartree(:)=dijhartree(:)-dij_ep(:)
         LIBPAW_DEALLOCATE(dij_ep)
       end if
       if (dijhartree_need) paw_ij(iatom)%dijhartree(:)=dijhartree(:)
     end if
     if (dij_need) then
       do idij=1,min(nsploop,2)
         if (cplex==1) then
           klmn1=1
           do klmn=1,lmn2_size
             paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+dijhartree(klmn)
             klmn1=klmn1+cplex_dij
           end do
         else
           paw_ij(iatom)%dij(:,idij)=paw_ij(iatom)%dij(:,idij)+dijhartree(:)
         end if
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
       LIBPAW_ALLOCATE(dijxc,(cplex_dij*lmn2_size,ndij))
       if (pawxcdev/=0) then
         LIBPAW_ALLOCATE(lmselect,(lm_size))
         lmselect(:)=paw_an(iatom)%lmselect(:)
         if (ipositron/=0) lmselect(:)=(lmselect(:).or.electronpositron_lmselect(:,iatom))
         call pawdijxcm(cplex,cplex_dij,dijxc,lmselect,ndij,nspden,nsppol,pawang,&
&                       pawrad(itypat),pawtab(itypat),paw_an(iatom)%vxc1,&
&                       paw_an(iatom)%vxct1,usexcnhat)
         LIBPAW_DEALLOCATE(lmselect)
       else
         call pawdijxc(cplex,cplex_dij,dijxc,ndij,nspden,nsppol,&
&                      pawang,pawrad(itypat),pawtab(itypat),paw_an(iatom)%vxc1,&
&                      paw_an(iatom)%vxct1,usexcnhat)
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
       LIBPAW_ALLOCATE(dijhat,(cplex_dij*lmn2_size,ndij))
       call pawdijhat(cplex,cplex_dij,dijhat,gprimd,iatom_tot,ipert,&
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
       if (dij_need) paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)+paw_ij(iatom)%dijnd(:,:)
     else

!    ===== Need to compute Dijnd
       LIBPAW_ALLOCATE(dijnd,(cplex_dij*lmn2_size,ndij))
       call pawdijnd(cplex_dij,dijnd,ndij,nucdipmom(:,iatom),pawrad(itypat),pawtab(itypat))
       if (dijnd_need) paw_ij(iatom)%dijnd(:,:)=dijnd(:,:)
       if (dij_need) paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)+dijnd(:,:)
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
       LIBPAW_ALLOCATE(dijso,(cplex_dij*lmn2_size,ndij))
       call pawdijso(cplex_dij,dijso,ndij,nspden,&
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
       lpawu=pawtab(itypat)%lpawu
       LIBPAW_ALLOCATE(dijpawu,(cplex_dij*lmn2_size,ndij))
       LIBPAW_POINTER_ALLOCATE(vpawu,(cplex_dij,lpawu*2+1,lpawu*2+1,nspden))
       if (pawtab(itypat)%usepawu>=10) vpawu=zero ! if dmft, do not apply U in LDA+U
       if (pawtab(itypat)%usepawu< 10) then
         call pawpupot(cplex_dij,ndij,paw_ij(iatom)%noccmmp,paw_ij(iatom)%nocctot,&
&                      nspden,pawprtvol,pawtab(itypat),vpawu)
       end if  
       if (natvshift_==0) then
         call pawdiju(cplex_dij,dijpawu,ndij,nspden,nsppol,pawtab(itypat),vpawu)
       else
         call pawdiju(cplex_dij,dijpawu,ndij,nspden,nsppol,pawtab(itypat),vpawu,&
&                     natvshift=natvshift_,atvshift=atvshift(:,:,iatom_tot),&
&                     fatvshift=fatvshift)
       end if
       LIBPAW_POINTER_DEALLOCATE(vpawu)
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
       if (dij_need) paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)+paw_ij(iatom)%dijexxc(:,:)
     else

!    ===== Need to compute DijEXXC
       LIBPAW_ALLOCATE(dijexxc,(cplex_dij*lmn2_size,ndij))
       if (pawxcdev/=0) then
         if (paw_ij(iatom)%has_exexch_pot/=2) then
           LIBPAW_POINTER_ALLOCATE(vpawx,(1,lmn2_size,nspden))
           call pawxpot(nspden,pawprtvol,pawrhoij(iatom),pawtab(itypat),vpawx)
         else
           vpawx=>paw_ij(iatom)%vpawx
         end if
         LIBPAW_ALLOCATE(lmselect,(lm_size))
         lmselect(:)=paw_an(iatom)%lmselect(:)
         if (ipositron/=0) lmselect(:)=(lmselect(:).or.electronpositron_lmselect(:,iatom))
         call pawdijexxc(cplex,cplex_dij,dijexxc,lmselect,ndij,nspden,nsppol,&
&             pawang,pawrad(itypat),pawtab(itypat),vpawx,paw_an(iatom)%vxc_ex)
         LIBPAW_DEALLOCATE(lmselect)
         if (paw_ij(iatom)%has_exexch_pot/=2) then
            LIBPAW_POINTER_DEALLOCATE(vpawx)
         end if
         if (dijexxc_need) paw_ij(iatom)%dijexxc(:,:)=dijexxc(:,:)
         if (dij_need) paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)+dijexxc(:,:)
         LIBPAW_DEALLOCATE(dijexxc)
       end if
     end if
 
   end if

!  ------------------------------------------------------------------------
!  ----------- Add Dij background contribution to the total Dij
!  ------------------------------------------------------------------------

   if (dij_need .AND. pawtab(itypat)%usepotzero==1 ) then
     do idij=1,min(nsploop,2)
       klmn1=1
       do klmn=1,lmn2_size
         paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+pawtab(itypat)%gammaij(klmn)*charge/ucvol
         klmn1=klmn1+cplex_dij
       end do
     end do
   end if


!  ------------------------------------------------------------------------
!  ----------- Compute Dijxc_hat
!  ------------------------------------------------------------------------

   if (dijxchat_need) then

     if (usexcnhat/=0) then
       LIBPAW_ALLOCATE(dijxchat,(cplex_dij*lmn2_size,ndij))
       call pawdijhat(cplex,cplex_dij,dijxchat,gprimd,iatom_tot,ipert,&
&                     natom,ndij,nfft,nfftot,nspden,nsppol,pawang,pawfgrtab(iatom),&
&                     pawtab(itypat),vxc,qphon,ucvol,xred,mpi_comm_grid=my_comm_grid)
       paw_ij(iatom)%dijxc_hat(:,:)=dijxchat(:,:)
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
       if (ipositron/=0) lmselect(:)=(lmselect(:).or.electronpositron_lmselect(:,iatom))
       call pawdijxcm(cplex,cplex_dij,dijxcval,lmselect,ndij,nspden,nsppol,&
&                     pawang,pawrad(itypat),pawtab(itypat),paw_an(iatom)%vxc1_val,&
&                     paw_an(iatom)%vxct1_val,0)
       LIBPAW_DEALLOCATE(lmselect)
     else
       call pawdijxc(cplex,cplex_dij,dijxcval,ndij,nspden,nsppol,&
&                    pawang,pawrad(itypat),pawtab(itypat),paw_an(iatom)%vxc1_val,&
&                    paw_an(iatom)%vxct1_val,0)
     end if
     paw_ij(iatom)%dijxc_val(:,:)=dijxcval(:,:)
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
!!  cplex=(RF calculations only) - 1 if RF 1st-order quantities are REAL, 2 if COMPLEX
!!  nspden=number of spin density components
!!  pawrhoij <type(pawrhoij_type)>= paw rhoij occupancies (and related data) for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data, for current atom
!!
!! OUTPUT
!!  dijxc(cplex*lmn2_size)=  D_ij^Hartree terms
!!
!! PARENTS
!!      m_pawdij,pawdenpot,pawdfptenergy
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijhartree(cplex,dijhartree,nspden,pawrhoij,pawtab)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdijhartree'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,nspden
!arrays
 real(dp),intent(out) :: dijhartree(:)
 type(pawrhoij_type),intent(in) :: pawrhoij
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: irhoij,ispden,jrhoij,kklmn,kklmn1,klmn,klmn1,lmn2_size,nspdiag
 character(len=500) :: msg
!arrays
 real(dp) :: ro(cplex)

! *************************************************************************

!Useful data
 lmn2_size=pawtab%lmn2_size
 nspdiag=1;if (nspden==2) nspdiag=2

!Check data consistency
 if (size(dijhartree,1)/=cplex*lmn2_size) then
   msg='invalid size for DijHartree !'
   MSG_BUG(msg)
 end if
 if (pawrhoij%cplex<cplex) then
   msg='  pawrhoij%cplex must be >=cplex  !'
   MSG_BUG(msg)
 end if

!------------------------------------------------------------------------
!----------- Allocations and initializations
!------------------------------------------------------------------------

 dijhartree=zero
!Real on-site quantities (ground-state calculation)
 if (cplex==1) then
   do ispden=1,nspdiag
     jrhoij=1
     do irhoij=1,pawrhoij%nrhoijsel
       klmn=pawrhoij%rhoijselect(irhoij)
       ro(1)=pawrhoij%rhoijp(jrhoij,ispden)*pawtab%dltij(klmn)
       dijhartree(klmn)=dijhartree(klmn)&
&       +ro(1)*pawtab%eijkl(klmn,klmn)
       do klmn1=1,klmn-1
         dijhartree(klmn1)=dijhartree(klmn1)&
&         +ro(1)*pawtab%eijkl(klmn1,klmn)
       end do
       do klmn1=klmn+1,lmn2_size
         dijhartree(klmn1)=dijhartree(klmn1)&
&         +ro(1)*pawtab%eijkl(klmn,klmn1)
       end do
       jrhoij=jrhoij+pawrhoij%cplex
     end do
   end do

!  Complex on-site quantities (response function calculation)
 else
   do ispden=1,nspdiag
     jrhoij=1
     do irhoij=1,pawrhoij%nrhoijsel
       klmn=pawrhoij%rhoijselect(irhoij);kklmn=2*klmn-1
       ro(1:2)=pawrhoij%rhoijp(jrhoij:jrhoij+1,ispden)*pawtab%dltij(klmn)
       dijhartree(kklmn:kklmn+1)=dijhartree(kklmn:kklmn+1)&
&       +ro(1:2)*pawtab%eijkl(klmn,klmn)
       do klmn1=1,klmn-1
         kklmn1=2*klmn1-1
         dijhartree(kklmn1:kklmn1+1)=dijhartree(kklmn1:kklmn1+1)&
&         +ro(1:2)*pawtab%eijkl(klmn1,klmn)
       end do
       do klmn1=klmn+1,lmn2_size
         kklmn1=2*klmn1-1
         dijhartree(kklmn1:kklmn1+1)=dijhartree(kklmn1:kklmn1+1)&
&         +ro(1:2)*pawtab%eijkl(klmn,klmn1)
       end do
       jrhoij=jrhoij+pawrhoij%cplex
     end do
   end do
 end if

end subroutine pawdijhartree
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
!!  cplex=(RF calculations only) - 1 if RF 1st-order quantities are REAL, 2 if COMPLEX
!!  cplex_dij=1 if dij is REAL, 2 if complex (2 for spin-orbit)
!!  ndij= number of spin components for Dij^hat
!!  nspden=number of spin density components
!!  nsppol=number of independent spin WF components
!!  pawang <type(pawang_type)>=paw angular mesh and related data, for current atom
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data, for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  vxc1(mesh_size,angl_size,nspden)=all-electron on-site XC potential for current atom
!!                                   given on a (r,theta,phi) grid
!!  vxct1(mesh_size,angl_size,nspden)=all-electron on-site XC potential for current atom
!!                                    given on a (r,theta,phi) grid
!!  usexcnhat= 1 if compensation density is included in Vxc, 0 otherwise
!!
!! OUTPUT
!!  dijxc(cplex_dij*lmn2_size,ndij)=  D_ij^XC terms
!!
!! NOTES
!!  cplex is for RF, cplex_dij is for non-collinear (nspinor==2)
!!
!! PARENTS
!!      m_pawdij
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijxc(cplex,cplex_dij,dijxc,ndij,nspden,nsppol,&
&                   pawang,pawrad,pawtab,vxc1,vxct1,usexcnhat)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdijxc'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,cplex_dij,ndij,nspden,nsppol,usexcnhat
 type(pawang_type),intent(in) :: pawang
!arrays
 real(dp),intent(in) :: vxc1(:,:,:),vxct1(:,:,:)
 real(dp),intent(out) :: dijxc(:,:)
 type(pawrad_type),intent(in) :: pawrad
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: angl_size,idij,idijend,ij_size,ilm,ils,ils1,ilslm,ipts,ir,ir1,isel,ispden
 integer :: jlm,j0lm,klmn,klmn1,klm,kln,l_size,lm0,lmax,lmin,lm_size,lmn2_size
 integer :: mesh_size,mm,nsploop
 real(dp) :: tmp,vi,vr,vxcijhat,vxcijhat_i
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: dijxc_idij(:),ff(:),gg(:),vxcij1(:),vxcij2(:),yylmr(:,:)

! *************************************************************************

!Useful data
 lm_size=pawtab%lcut_size**2
 lmn2_size=pawtab%lmn2_size
 ij_size=pawtab%ij_size
 l_size=pawtab%l_size
 mesh_size=pawtab%mesh_size
 angl_size=pawang%angl_size

!Check data consistency
 if (size(dijxc,1)/=cplex_dij*lmn2_size.or.size(dijxc,2)/=ndij) then
   msg='invalid sizes for Dijxc !'
   MSG_BUG(msg)
 end if
 if (size(vxc1,1)/=cplex*mesh_size.or.size(vxct1,1)/=cplex*mesh_size.or.&
&    size(vxc1,2)/=angl_size.or.size(vxct1,2)/=angl_size.or.&
&    size(vxc1,3)/=nspden.or.size(vxct1,3)/=nspden) then
   msg='invalid sizes for vxc1 or vxct1 !'
   MSG_BUG(msg)
 end if
 if (cplex_dij<cplex) then
   msg='cplex_dij must be >= cplex !'
   MSG_BUG(msg)
 end if

!Precompute products Ylm*Ylpmp
 lmax=maxval(pawtab%indklmn(4,1:lmn2_size))
 LIBPAW_ALLOCATE(yylmr,(lmax**2*(lmax**2+1)/2,angl_size))
 do ipts=1,angl_size
   do jlm=1,lmax**2
     j0lm=jlm*(jlm-1)/2
     do ilm=1,jlm
       klm=j0lm+ilm
       yylmr(klm,ipts)=pawang%ylmr(ilm,ipts)*pawang%ylmr(jlm,ipts)
     end do
   end do
 end do

!Init memory
 dijxc=zero
 LIBPAW_ALLOCATE(dijxc_idij,(cplex*lmn2_size))
 LIBPAW_ALLOCATE(vxcij1,(cplex*ij_size))
 LIBPAW_ALLOCATE(vxcij2,(cplex*l_size))
 LIBPAW_ALLOCATE(ff,(mesh_size))
 LIBPAW_ALLOCATE(gg,(mesh_size))

!----------------------------------------------------------
!Loop over spin components
!----------------------------------------------------------
 nsploop=nsppol;if (ndij==4) nsploop=4
 do idij=1,nsploop
   if (idij<=nsppol.or.(nspden==4.and.idij<=3).or.cplex==2) then

     idijend=idij+idij/3;if (cplex==2) idijend=idij
     do ispden=idij,idijend

       dijxc_idij=zero

!      ----------------------------------------------------------
!      Loop on angular mesh
!      ----------------------------------------------------------
       do ipts=1,angl_size

!        ===== Vxc_ij_1 (tmp) =====
         vxcij1=zero
         if (cplex==1) then
           do kln=1,ij_size
             ff(1:mesh_size)= &
&               vxc1(1:mesh_size,ipts,ispden)*pawtab%phiphj(1:mesh_size,kln) &
&              -vxct1(1:mesh_size,ipts,ispden)*pawtab%tphitphj(1:mesh_size,kln)
             call simp_gen(vxcij1(kln),ff,pawrad)
           end do
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
         end if

!        ===== Vxc_ij_2 (tmp) =====
         vxcij2=zero
         if (usexcnhat/=0) then
           if (cplex==1) then
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
         if (cplex==1) then
           do klmn=1,lmn2_size
             klm=pawtab%indklmn(1,klmn);kln=pawtab%indklmn(2,klmn)
             lmin=pawtab%indklmn(3,klmn);lmax=pawtab%indklmn(4,klmn)
             dijxc_idij(klmn)=dijxc_idij(klmn) &
&                            +vxcij1(kln)*pawang%angwgth(ipts)*yylmr(klm,ipts)*four_pi
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
             klmn1=klmn1+cplex
           end do ! Loop klmn
         end if

!      ----------------------------------------------------------
!      End loop on angular points
       end do

!      ----------------------------------------------------------
!      Deduce some part of Dij according to symmetries
!      ----------------------------------------------------------

       if (cplex==1) then
         if (ispden<3) then
           if (cplex_dij==1) then
             dijxc(1:lmn2_size,idij)=dijxc_idij(1:lmn2_size)
           else
             klmn1=1
             do klmn=1,lmn2_size
               dijxc(klmn1  ,idij)=dijxc_idij(klmn)
               dijxc(klmn1+1,idij)=zero
               klmn1=klmn1+cplex_dij
             end do
           end if
         else
           klmn1=max(1,ispden-2)
           do klmn=1,lmn2_size
             dijxc(klmn1,idij)=dijxc_idij(klmn)
             klmn1=klmn1+cplex_dij
           end do
         end if
       else !cplex=2
         if (ispden<=3) then
           dijxc(1:cplex*lmn2_size,idij)=dijxc_idij(1:cplex*lmn2_size)
         else     
           klmn1=1  ! Remember V(4) contains i.V^21
           do klmn=1,lmn2_size
             dijxc(klmn1  ,idij)= dijxc_idij(klmn+1)
             dijxc(klmn1+1,idij)=-dijxc_idij(klmn  )
             klmn1=klmn1+cplex_dij
           end do
         end if
       end if

     end do !ispden

   else if (nspden==4.and.idij==4) then ! cplex=1 here
     dijxc(:,idij)=dijxc(:,idij-1)
     klmn1=2
     do klmn=1,lmn2_size
       dijxc(klmn1,idij)=-dijxc(klmn1,idij)
       klmn1=klmn1+cplex_dij
     end do
   else if (nsppol==1.and.idij==2) then ! cplex=1 here
     dijxc(:,idij)=dijxc(:,idij-1)
   end if

!----------------------------------------------------------
!End loop on spin density components
 end do

!Free temporary memory spaces
 LIBPAW_DEALLOCATE(yylmr)
 LIBPAW_DEALLOCATE(dijxc_idij)
 LIBPAW_DEALLOCATE(vxcij1)
 LIBPAW_DEALLOCATE(vxcij2)
 LIBPAW_DEALLOCATE(ff)
 LIBPAW_DEALLOCATE(gg)

end subroutine pawdijxc
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
!!  cplex=(RF calculations only) - 1 if RF 1st-order quantities are REAL, 2 if COMPLEX
!!  cplex_dij=1 if dij is REAL, 2 if complex (2 for spin-orbit)
!!  [hyb_mixing, hyb_mixing_sr]= -- optional-- mixing factors for the global (resp. screened) XC hybrid functional
!!  ndij= number of spin components for Dij^hat
!!  nspden=number of spin density components
!!  pawrhoij <type(pawrhoij_type)>= paw rhoij occupancies (and related data) for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data, for current atom
!!  cplex=(RF calculations only) - 1 if RF 1st-order quantities are REAL, 2 if COMPLEX
!!
!!  lmselect(lm_size)=select the non-zero LM-moments of on-site potentials
!!  nspden=number of spin density components
!!  nsppol=number of independent spin WF components
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data, for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data, for current atom
!!
!! OUTPUT
!!  dijfock_vv(cplex_dij*lmn2_size,ndij)=  D_ij^fock terms for valence-valence interactions
!!  dijfock_cv(cplex_dij*lmn2_size,ndij)=  D_ij^fock terms for core-valence interactions
!!
!! PARENTS
!!      m_pawdij,pawdenpot
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijfock(cplex,cplex_dij,dijfock_vv,dijfock_cv,hyb_mixing,hyb_mixing_sr,ndij,nspden,nsppol,pawrhoij,pawtab)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdijfock'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,cplex_dij,ndij,nspden,nsppol
 real(dp),intent(in) :: hyb_mixing,hyb_mixing_sr
!arrays
 real(dp),intent(out) :: dijfock_vv(:,:),dijfock_cv(:,:)
 type(pawrhoij_type),intent(in) :: pawrhoij
 type(pawtab_type),intent(in),target :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: idij,idijend,ispden,irhoij,jrhoij,ilmn_i,jlmn_j,ilmn_k,jlmn_l
 integer :: klmn_kl,klmn_ij,klmn_il,klmn_kj,klmn,klmn1,nsploop,lmn2_size
 
 character(len=500) :: msg
!arrays
 real(dp) :: ro(cplex)
 real(dp),allocatable :: dijfock_idij_vv(:),dijfock_idij_cv(:)
 real(dp),pointer :: eijkl(:,:)
  
! *************************************************************************

!Useful data
 lmn2_size=pawtab%lmn2_size

!Check data consistency
 if (size(dijfock_vv,1)/=cplex_dij*lmn2_size.or.size(dijfock_vv,2)/=ndij) then
   msg='invalid sizes for Dijfock_vv !'
   MSG_BUG(msg)
 end if
 if (size(dijfock_cv,1)/=cplex_dij*lmn2_size.or.size(dijfock_cv,2)/=ndij) then
   msg='invalid sizes for Dijfock_cv !'
   MSG_BUG(msg)
 end if
 if (cplex_dij<cplex) then
   msg='cplex_dij must be >= cplex !'
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
 LIBPAW_ALLOCATE(dijfock_idij_vv,(cplex*lmn2_size))
 LIBPAW_ALLOCATE(dijfock_idij_cv,(cplex*lmn2_size))

!----------------------------------------------------------
!Loop over spin components
!----------------------------------------------------------
 nsploop=nsppol;if (ndij==4) nsploop=4
 do idij=1,nsploop
   if (idij<=nsppol.or.(nspden==4.and.idij<=3).or.cplex==2) then

     idijend=idij+idij/3;if (cplex==2) idijend=idij
     do ispden=idij,idijend
!!!! WARNING : What follows has been tested only for cases where nsppol=1 and 2, nspden=1 and 2 with nspinor=1. 
       dijfock_idij_vv=zero
       dijfock_idij_cv=zero
!Real on-site quantities (ground-state calculation)
       if (cplex==1) then
!* Loop on the non-zero elements rho_kl
         do irhoij=1,pawrhoij%nrhoijsel
           klmn_kl=pawrhoij%rhoijselect(irhoij)
           ro(1)=pawrhoij%rhoijp(irhoij,ispden)*pawtab%dltij(klmn_kl)
           ilmn_k=pawtab%indklmn(7,klmn_kl)
           jlmn_l=pawtab%indklmn(8,klmn_kl)
!* Fock contribution to the element (k,l) of dijfock 
           dijfock_idij_vv(klmn_kl)=dijfock_idij_vv(klmn_kl)-ro(1)*eijkl(klmn_kl,klmn_kl)
!* Fock contribution to the element (i,j) of dijfock with (i,j) < (k,l)
!* We remind that i<j and k<l by construction
           do klmn_ij=1,klmn_kl-1
             ilmn_i=pawtab%indklmn(7,klmn_ij)
             jlmn_j=pawtab%indklmn(8,klmn_ij)
!* In this case, i < l
             klmn_il=jlmn_l*(jlmn_l-1)/2+ilmn_i
!* If k >j, one must consider the index of the symmetric element (j,k) ; otherwise, the index of the element (k,j) is calculated. 
             if (ilmn_k>jlmn_j) then
               klmn_kj=ilmn_k*(ilmn_k-1)/2+jlmn_j
             else
               klmn_kj=jlmn_j*(jlmn_j-1)/2+ilmn_k
             end if
!* In this case, (i,l) >= (k,j)
             dijfock_idij_vv(klmn_ij)=dijfock_idij_vv(klmn_ij)-ro(1)*eijkl(klmn_il,klmn_kj)
           end do
!* Fock contribution to the element (i,j) of dijfock with (i,j) > (k,l)
!* We remind that i<j and k<l by construction           
           do klmn_ij=klmn_kl+1,lmn2_size
             ilmn_i=pawtab%indklmn(7,klmn_ij)
             jlmn_j=pawtab%indklmn(8,klmn_ij)
!* In this case, k < j
             klmn_kj=jlmn_j*(jlmn_j-1)/2+ilmn_k
!* If i >l, one must consider the index of the symmetric element (l,i) ; otherwise, the index of the element (i,l) is calculated. 
             if (ilmn_i>jlmn_l) then
               klmn_il=ilmn_i*(ilmn_i-1)/2+jlmn_l
             else
               klmn_il=jlmn_l*(jlmn_l-1)/2+ilmn_i
             end if
!* In this case, (k,j) >= (i,l)
             dijfock_idij_vv(klmn_ij)=dijfock_idij_vv(klmn_ij)-ro(1)*eijkl(klmn_kj,klmn_il)
           end do
         end do
! Add the core-valence contribution
         do klmn_ij=1,lmn2_size
           dijfock_idij_cv(klmn_ij)=dijfock_idij_cv(klmn_ij)+pawtab%ex_cvij(klmn_ij)
         end do  

!Complex on-site quantities (response function calculation)
       else !cplex=2
         jrhoij=1
!* Loop on the non-zero elements rho_kl         
         do irhoij=1,pawrhoij%nrhoijsel
           klmn_kl=pawrhoij%rhoijselect(irhoij)
           ro(1)=pawrhoij%rhoijp(jrhoij,ispden)*pawtab%dltij(klmn_kl)
           ro(2)=pawrhoij%rhoijp(jrhoij+1,ispden)*pawtab%dltij(klmn_kl)
           ilmn_k=pawtab%indklmn(7,klmn_kl)
           jlmn_l=pawtab%indklmn(8,klmn_kl)
!* Fock contribution to the element (k,l) of dijfock            
           dijfock_idij_vv(klmn_kl)=dijfock_idij_vv(klmn_kl)-ro(1)*eijkl(klmn_kl,klmn_kl)
           dijfock_idij_vv(klmn_kl+1)=dijfock_idij_vv(klmn_kl)-ro(2)*eijkl(klmn_kl,klmn_kl)
!* Fock contribution to the element (i,j) of dijfock with (i,j) < (k,l)
!* We remind that i<j and k<l by construction
           do klmn_ij=1,klmn_kl-1
             ilmn_i=pawtab%indklmn(7,klmn_ij)
             jlmn_j=pawtab%indklmn(8,klmn_ij)
!* In this case, i < l
             klmn_il=jlmn_l*(jlmn_l-1)/2+ilmn_i
!* If k >j, one must consider the index of the symmetric element (j,k) ; otherwise, the index of the element (k,j) is calculated. 
             if (ilmn_k>jlmn_j) then
               klmn_kj=ilmn_k*(ilmn_k-1)/2+jlmn_j
             else
               klmn_kj=jlmn_j*(jlmn_j-1)/2+ilmn_k
             end if
!* In this case, (i,l) >= (k,j)
             dijfock_idij_vv(klmn_ij)=dijfock_idij_vv(klmn_ij)-ro(1)*eijkl(klmn_il,klmn_kj)
             dijfock_idij_vv(klmn_ij+1)=dijfock_idij_vv(klmn_ij)-ro(2)*eijkl(klmn_il,klmn_kj)
           end do
!* Fock contribution to the element (i,j) of dijfock with (i,j) > (k,l)
!* We remind that i<j and k<l by construction           
           do klmn_ij=klmn_kl+1,lmn2_size
             ilmn_i=pawtab%indklmn(7,klmn_ij)
             jlmn_j=pawtab%indklmn(8,klmn_ij)
!* In this case, k < j
             klmn_kj=jlmn_j*(jlmn_j-1)/2+ilmn_k
!* If i >l, one must consider the index of the symmetric element (l,i) ; otherwise, the index of the element (i,l) is calculated. 
             if (ilmn_i>jlmn_l) then
               klmn_kj=ilmn_i*(ilmn_i-1)/2+jlmn_l
             else
               klmn_kj=jlmn_l*(jlmn_l-1)/2+ilmn_i
             end if
!* In this case, (k,j) >= (i,l)
             dijfock_idij_vv(klmn_ij)=dijfock_idij_vv(klmn_ij)-ro(1)*eijkl(klmn_kj,klmn_il)
             dijfock_idij_vv(klmn_ij+1)=dijfock_idij_vv(klmn_ij)-ro(2)*eijkl(klmn_kj,klmn_il)
           end do

           jrhoij=jrhoij+cplex
         end do 
! Add the core-valence contribution
         do klmn_ij=1,lmn2_size,2
           dijfock_idij_cv(klmn_ij)=dijfock_idij_cv(klmn_ij)+pawtab%ex_cvij(klmn_ij)
         end do  
         
       end if

!      ----------------------------------------------------------
!      Deduce some part of Dij according to symmetries
!      ----------------------------------------------------------

       if (cplex==1) then
         if (ispden<3) then
           if (cplex_dij==1) then
             dijfock_vv(1:lmn2_size,idij)=dijfock_idij_vv(1:lmn2_size)
             dijfock_cv(1:lmn2_size,idij)=dijfock_idij_cv(1:lmn2_size)
           else
             klmn1=1
             do klmn=1,lmn2_size
               dijfock_vv(klmn1  ,idij)=dijfock_idij_vv(klmn)
               dijfock_cv(klmn1  ,idij)=dijfock_idij_cv(klmn)
               dijfock_vv(klmn1+1,idij)=zero
               dijfock_cv(klmn1+1,idij)=zero
               klmn1=klmn1+cplex_dij
             end do
           end if
         else
           klmn1=max(1,ispden-2)
           do klmn=1,lmn2_size
             dijfock_vv(klmn1,idij)=dijfock_idij_vv(klmn)
             dijfock_cv(klmn1,idij)=dijfock_idij_cv(klmn)
             klmn1=klmn1+cplex_dij
           end do
         end if
       else !cplex=2
         if (ispden<=3) then
           dijfock_vv(1:cplex*lmn2_size,idij)=dijfock_idij_vv(1:cplex*lmn2_size)
           dijfock_cv(1:cplex*lmn2_size,idij)=dijfock_idij_cv(1:cplex*lmn2_size)
         else     
           klmn1=1  ! Remember V(4) contains i.V^21
           do klmn=1,lmn2_size
             dijfock_vv(klmn1  ,idij)= dijfock_idij_vv(klmn+1)
             dijfock_vv(klmn1+1,idij)=-dijfock_idij_vv(klmn  )
             dijfock_cv(klmn1  ,idij)= dijfock_idij_cv(klmn+1)
             dijfock_cv(klmn1+1,idij)=-dijfock_idij_cv(klmn  )
             klmn1=klmn1+cplex_dij
           end do
         end if
       end if

     end do !ispden

   else if (nspden==4.and.idij==4) then ! cplex=1 here
     dijfock_vv(:,idij)=dijfock_vv(:,idij-1)
     dijfock_cv(:,idij)=dijfock_cv(:,idij-1)
     klmn1=2
     do klmn=1,lmn2_size
       dijfock_vv(klmn1,idij)=-dijfock_vv(klmn1,idij)
       dijfock_cv(klmn1,idij)=-dijfock_cv(klmn1,idij)
       klmn1=klmn1+cplex_dij
     end do
   else if (nsppol==1.and.idij==2) then ! cplex=1 here
     dijfock_vv(:,idij)=dijfock_vv(:,idij-1)
     dijfock_cv(:,idij)=dijfock_cv(:,idij-1)
   end if

!----------------------------------------------------------
!End loop on spin density components
 end do

 if (abs(hyb_mixing)>tol8) then
   dijfock_vv(:,:) = hyb_mixing*dijfock_vv(:,:)
 else if (abs(hyb_mixing_sr)>tol8) then
   dijfock_vv(:,:) = hyb_mixing_sr*dijfock_vv(:,:)
 end if
 dijfock_cv(:,:) = (hyb_mixing+hyb_mixing_sr)*dijfock_cv(:,:)

!Free temporary memory spaces
 LIBPAW_DEALLOCATE(dijfock_idij_vv)
 LIBPAW_DEALLOCATE(dijfock_idij_cv)

end subroutine pawdijfock
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
!!  cplex=(RF calculations only) - 1 if RF 1st-order quantities are REAL, 2 if COMPLEX
!!  cplex_dij=1 if dij is REAL, 2 if complex (2 for spin-orbit)
!!  lmselect(lm_size)=select the non-zero LM-moments of on-site potentials
!!  ndij= number of spin components for Dij^hat
!!  nspden=number of spin density components
!!  nsppol=number of independent spin WF components
!!  pawang <type(pawang_type)>=paw angular mesh and related data, for current atom
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data, for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  vxc1(mesh_size,lm_size,nspden)=all-electron on-site XC potential for current atom
!!                                 given on (l,m) spherical moments
!!  vxct1(mesh_size,lm_size,nspden)=all-electron on-site XC potential for current atom
!!                                  given on (l,m) spherical moments
!!  usexcnhat= 1 if compensation density is included in Vxc, 0 otherwise
!!
!! OUTPUT
!!  dijxc(cplex_dij*lmn2_size,ndij)=  D_ij^XC terms
!!
!! NOTES
!!  cplex is for RF, cplex_dij is for non-collinear (nspinor==2)
!!
!! PARENTS
!!      m_pawdij
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijxcm(cplex,cplex_dij,dijxc,lmselect,ndij,nspden,nsppol,&
&                    pawang,pawrad,pawtab,vxc1,vxct1,usexcnhat)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdijxcm'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,cplex_dij,ndij,nspden,nsppol,usexcnhat
 type(pawang_type),intent(in) :: pawang
!arrays
 logical :: lmselect(:)
 real(dp),intent(in) :: vxc1(:,:,:),vxct1(:,:,:)
 real(dp),intent(out) :: dijxc(:,:)
 type(pawrad_type),intent(in) :: pawrad
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: idij,idijend,ij_size,ir,ir1,isel,ispden,klm,klm1,klmn,klmn1,kln
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
 if (size(dijxc,1)/=cplex_dij*lmn2_size.or.size(dijxc,2)/=ndij) then
   msg='invalid sizes for Dijxc !'
   MSG_BUG(msg)
 end if
 if (size(lmselect)/=lm_size) then
   msg='invalid size for lmselect !'
   MSG_BUG(msg)
 end if
 if (size(vxc1,1)/=cplex*mesh_size.or.size(vxct1,1)/=cplex*mesh_size.or.&
&    size(vxc1,2)/=lm_size.or.size(vxct1,2)/=lm_size.or.&
&    size(vxc1,3)/=nspden.or.size(vxct1,3)/=nspden) then
   msg='invalid sizes for vxc1 or vxct1 !'
   MSG_BUG(msg)
 end if
 if (cplex_dij<cplex) then
   msg='cplex_dij must be >= cplex !'
   MSG_BUG(msg)
 end if

!Init memory
 dijxc=zero
 LIBPAW_ALLOCATE(dijxc_idij,(cplex*lmn2_size))
 LIBPAW_ALLOCATE(vxcij1,(cplex*ij_size))
 LIBPAW_ALLOCATE(ff,(mesh_size))
 LIBPAW_ALLOCATE(gg,(mesh_size))

!----------------------------------------------------------
!Loop over spin components
!----------------------------------------------------------
 nsploop=nsppol;if (ndij==4) nsploop=4
 do idij=1,nsploop
   if (idij<=nsppol.or.(nspden==4.and.idij<=3).or.cplex==2) then

     idijend=idij+idij/3;if (cplex==2) idijend=idij
     do ispden=idij,idijend

       dijxc_idij=zero

!      ----------------------------------------------------------
!      Summing over (l,m) moments
!      ----------------------------------------------------------
       do klm=1,lm_size
         if (lmselect(klm)) then

!          ===== Vxc_ij_1 (tmp) =====
           vxcij1=zero
           if (cplex==1) then
             do kln=1,ij_size
               ff(1:mesh_size)= &
&                 vxc1(1:mesh_size,klm,ispden)*pawtab%phiphj(1:mesh_size,kln) &
&                -vxct1(1:mesh_size,klm,ispden)*pawtab%tphitphj(1:mesh_size,kln)
               call simp_gen(vxcij1(kln),ff,pawrad)
             end do
           else ! cplex==2
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
             if (cplex==1) then
               ff(1:mesh_size)=vxct1(1:mesh_size,klm,ispden) &
&                             *pawtab%shapefunc(1:mesh_size,ll) &
&                             *pawrad%rad(1:mesh_size)**2
               call simp_gen(vxcij2,ff,pawrad)
             else ! cplex==2
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
           if (cplex==1) then
             do klmn=1,lmn2_size
               klm1=pawtab%indklmn(1,klmn)
               kln=pawtab%indklmn(2,klmn)
               isel=pawang%gntselect(klm,klm1)
               if (isel>0) &
&                dijxc_idij(klmn)=dijxc_idij(klmn)+vxcij1(kln)*pawang%realgnt(isel)
               if (usexcnhat/=0) &
                 dijxc_idij(klmn)=dijxc_idij(klmn)-pawtab%qijl(klm,klmn)*vxcij2
             end do ! Loop klmn
           else ! cplex==2
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
               klmn1=klmn1+cplex
             end do ! Loop klmn
           end if

         end if ! klm selection
       end do  ! Loop klm

!      ----------------------------------------------------------
!      Deduce some part of Dij according to symmetries
!      ----------------------------------------------------------

       if (cplex==1) then
         if (ispden<3) then
           if (cplex_dij==1) then
             dijxc(1:lmn2_size,idij)=dijxc_idij(1:lmn2_size)
           else
             klmn1=1
             do klmn=1,lmn2_size
               dijxc(klmn1  ,idij)=dijxc_idij(klmn)
               dijxc(klmn1+1,idij)=zero
               klmn1=klmn1+cplex_dij
             end do
           end if
         else
           klmn1=max(1,ispden-2)
           do klmn=1,lmn2_size
             dijxc(klmn1,idij)=dijxc_idij(klmn)
             klmn1=klmn1+cplex_dij
           end do
         end if
       else !cplex=2
         if (ispden<=3) then
           dijxc(1:cplex*lmn2_size,idij)=dijxc_idij(1:cplex*lmn2_size)
         else     
           klmn1=1  ! Remember V(4) contains i.V^21
           do klmn=1,lmn2_size
             dijxc(klmn1  ,idij)= dijxc_idij(klmn+1)
             dijxc(klmn1+1,idij)=-dijxc_idij(klmn  )
             klmn1=klmn1+cplex_dij
           end do
         end if
       end if

     end do !ispden

   else if (nspden==4.and.idij==4) then ! cplex=1 here
     dijxc(:,idij)=dijxc(:,idij-1)
     klmn1=2
     do klmn=1,lmn2_size
       dijxc(klmn1,idij)=-dijxc(klmn1,idij)
       klmn1=klmn1+cplex_dij
     end do
   else if (nsppol==1.and.idij==2) then ! cplex=1 here
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
!!  cplex=(RF calculations only) - 1 if RF 1st-order quantities are REAL, 2 if COMPLEX
!!  cplex_dij=1 if dij is REAL, 2 if complex (2 for spin-orbit)
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!  iatom=absolute index of current atom (between 1 and natom)
!!  ipert=index of perturbation; used only for RF calculation ; set ipert<=0 for GS calculations.
!!  natom=total number of atoms
!!  ndij= number of spin components for Dij^hat
!!  ngrid=number of points of the real space grid (FFT, WVL, ...) treated by current proc
!!  ngridtot=total number of points of the real space grid (FFT, WVL, ...)
!!           For the FFT grid, thi should be equal to ngfft1*ngfft2*ngfft3
!!  nspden=number of spin density components
!!  nsppol=number of independent spin WF components
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab<type(pawfgrtab_type)>=atomic data given on fine rectangular grid for current atom
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data, for current atom
!!  Pot(cplex*ngrid,nspden)=potential on real space grid
!!  qphon(3)=(RF calculations only) - wavevector of the phonon
!!  ucvol=unit cell volume
!!  xred(3,my_natom)= reduced atomic coordinates
!!
!! OUTPUT
!!  dijhat(cplex_dij*lmn2_size,ndij)= D_ij^hat terms
!!
!! NOTES
!!  cplex is for RF, cplex_dij is for non-collinear (nspinor==2)
!!
!! PARENTS
!!      fock_getghc,m_pawdij
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijhat(cplex,cplex_dij,dijhat,gprimd,iatom,ipert,&
&                    natom,ndij,ngrid,ngridtot,nspden,nsppol,pawang,pawfgrtab,&
&                    pawtab,Pot,qphon,ucvol,xred,&
&                    mpi_comm_grid) ! Optional argument


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdijhat'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,cplex_dij,iatom,ipert,natom,ndij
 integer,intent(in) :: ngrid,ngridtot,nspden,nsppol
 integer,intent(in),optional :: mpi_comm_grid
 real(dp),intent(in) :: ucvol
 type(pawang_type),intent(in) :: pawang
 type(pawfgrtab_type),intent(inout) :: pawfgrtab
!arrays
 real(dp),intent(in) :: gprimd(3,3),Pot(cplex*ngrid,nspden),qphon(3),xred(3,natom)
 real(dp),intent(out) :: dijhat(:,:)
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: ic,idij,idijend,ier,ils,ilslm,ilslm1,isel,ispden,jc,klm,klmn,klmn1
 integer :: lm0,lm_size,lmax,lmin,lmn2_size,mm,my_comm_grid,nfgd,nsploop,optgr0,optgr1
 logical :: has_phase,qne0
 real(dp) :: vi,vr
 character(len=500) :: msg
!arrays
 real(dp) :: rdum(1)
 real(dp),allocatable :: dijhat_idij(:),prod(:)

! *************************************************************************

!Useful data
 lm_size=pawtab%lcut_size**2
 lmn2_size=pawtab%lmn2_size
 nfgd=pawfgrtab%nfgd
 has_phase=.false.
 qne0=(qphon(1)**2+qphon(2)**2+qphon(3)**2>=1.d-15)
 my_comm_grid=xmpi_comm_self;if (present(mpi_comm_grid)) my_comm_grid=mpi_comm_grid


!Check data consistency
 if (size(dijhat,1)/=cplex_dij*lmn2_size.or.size(dijhat,2)/=ndij) then
   msg='invalid sizes for Dijhat !'
   MSG_BUG(msg)
 end if
 if (pawfgrtab%rfgd_allocated==0.and.ipert>0.and.ipert<=natom.and.qne0) then
   msg='pawfgrtab()%rfgd array must be allocated  !'
   MSG_BUG(msg)
 end if
 if (cplex_dij<cplex) then
   msg='cplex_dij must be >= cplex !'
   MSG_BUG(msg)
 end if

!Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
 if ((pawfgrtab%gylm_allocated==0).or.((ipert==iatom).and.(pawfgrtab%gylmgr_allocated==0))) then
   optgr0=0;optgr1=0
   if (pawfgrtab%gylm_allocated==0) then
     if (allocated(pawfgrtab%gylm))  then
       LIBPAW_DEALLOCATE(pawfgrtab%gylm)
     end if
     LIBPAW_ALLOCATE(pawfgrtab%gylm,(nfgd,lm_size))
     pawfgrtab%gylm_allocated=2;optgr0=1
   end if
   if ((ipert==iatom).and.(pawfgrtab%gylmgr_allocated==0)) then
     if (allocated(pawfgrtab%gylmgr))  then
       LIBPAW_DEALLOCATE(pawfgrtab%gylmgr)
     end if
     LIBPAW_ALLOCATE(pawfgrtab%gylmgr,(3,nfgd,lm_size))
     pawfgrtab%gylmgr_allocated=2;optgr1=1
   end if
   if (optgr0+optgr1>0) then
     call pawgylm(pawfgrtab%gylm,pawfgrtab%gylmgr,rdum,lm_size,nfgd,optgr0,optgr1,0,&
&                 pawtab,pawfgrtab%rfgd)
   end if
 end if

!Eventually compute exp(i.q.r) factors for the current atom (if not already done)
 if ((ipert==iatom).and.qne0.and.(pawfgrtab%expiqr_allocated==0)) then
   if (allocated(pawfgrtab%expiqr))  then
     LIBPAW_DEALLOCATE(pawfgrtab%expiqr)
   end if
   LIBPAW_ALLOCATE(pawfgrtab%expiqr,(2,nfgd))
   call pawexpiqr(pawfgrtab%expiqr,gprimd,nfgd,qphon,pawfgrtab%rfgd,xred(:,iatom))
   pawfgrtab%expiqr_allocated=2
 end if
 has_phase=(qne0.and.ipert>0.and.pawfgrtab%expiqr_allocated/=0)

!Init memory
 dijhat=zero
 LIBPAW_ALLOCATE(prod,(cplex*lm_size))
 LIBPAW_ALLOCATE(dijhat_idij,(cplex*lmn2_size))

!----------------------------------------------------------
!Loop over spin components
!----------------------------------------------------------
 nsploop=nsppol;if (ndij==4) nsploop=4
 do idij=1,nsploop
   if (idij<=nsppol.or.(nspden==4.and.idij<=3).or.cplex==2) then

     idijend=idij+idij/3;if (cplex==2) idijend=idij
     do ispden=idij,idijend

!      ------------------------------------------------------
!      Compute Int[V(r).g_l(r).Y_lm(r)]
!      ------------------------------------------------------
!       Note for non-collinear magnetism:
!          We compute Int[V^(alpha,beta)(r).g_l(r).Y_lm(r)]
!          Remember: if nspden=4, V is stored as : V^11, V^22, V^12, i.V^21

       prod=zero

!      ===== Standard case ============================
       if (.not.has_phase) then
         if (cplex==1) then
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
             ilslm1=ilslm1+cplex
           end do
         end if
 
!      ===== Including Exp(iqr) phase (DFPT only) =====
       else
         if (cplex==1) then
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
             ilslm1=ilslm1+cplex
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

       if (cplex==1) then
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
 
       if (cplex==1) then
         if (ispden<3) then
           if (cplex_dij==1) then
             dijhat(1:lmn2_size,idij)=dijhat_idij(1:lmn2_size)
           else
             klmn1=1
             do klmn=1,lmn2_size
               dijhat(klmn1  ,idij)=dijhat_idij(klmn)
               dijhat(klmn1+1,idij)=zero
               klmn1=klmn1+cplex_dij
             end do
           end if
         else
           klmn1=max(1,ispden-2)
           do klmn=1,lmn2_size
             dijhat(klmn1,idij)=dijhat_idij(klmn)
             klmn1=klmn1+cplex_dij
           end do
         end if
       else !cplex=2
         if (ispden<=3) then
           dijhat(1:cplex*lmn2_size,idij)=dijhat_idij(1:cplex*lmn2_size)
         else     
           klmn1=1  ! Remember V(4) contains i.V^21
           do klmn=1,lmn2_size
             dijhat(klmn1  ,idij)= dijhat_idij(klmn+1)
             dijhat(klmn1+1,idij)=-dijhat_idij(klmn  )
             klmn1=klmn1+cplex_dij
           end do
         end if
       end if

     end do !ispden

   else if (nspden==4.and.idij==4) then ! cplex=1 here
     dijhat(:,idij)=dijhat(:,idij-1)
     klmn1=2
     do klmn=1,lmn2_size
       dijhat(klmn1,idij)=-dijhat(klmn1,idij)
       klmn1=klmn1+cplex_dij
     end do
   else if (nsppol==1.and.idij==2) then ! cplex=1 here
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
 if (pawfgrtab%gylmgr_allocated==2) then
   LIBPAW_DEALLOCATE(pawfgrtab%gylmgr)
   LIBPAW_ALLOCATE(pawfgrtab%gylmgr,(0,0,0))
   pawfgrtab%gylmgr_allocated=0
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
!!  ndij= number of spin components
!!  nucdipmom(3) nuclear magnetic dipole moment for current atom
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data, for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data, for current atom
!!
!! OUTPUT
!!  dijnd(cplex_dij*lmn2_size,ndij)= nuclear dipole moment Dij terms
!!  cplex_dij=2 must be 2
!!
!! NOTES
!!   On-site contribution of a nuclear magnetic dipole moment at $R$. Hamiltonian is
!!   $H=(1/2m_e)(p - q_e A)^2 + V$, and vector potential $A$ is
!!   $A=(\mu_0/4\pi) m\times (r-R)/|r-R|^3 = (\mu_0/4\pi) L_R\cdot m/|r-R|^3$ where
!!   $L_R$ is the on-site orbital angular momentum and $m$ is the nuclear magnetic 
!!   dipole moment. For an electron (as usual), mass m_e = 1 and charge q_e = -1.
!!   Second order term in A is ignored.
!!
!! PARENTS
!!      m_pawdij,pawdenpot
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijnd(cplex_dij,dijnd,ndij,nucdipmom,pawrad,pawtab)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdijnd'
!End of the abilint section

 implicit none

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
 integer :: idir,ilmn,il,im,iln,ilm,jlmn,jl,jm,jlm,jln,j0lmn,klmn,kln,mesh_size
 real(dp) :: intgr3,permeability
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

! magnetic permeability mu_0/four_pi in atomic units
! this constant is also used in getghcnd.F90, if you change it here,
! change it there also for consistency
 permeability=5.325135453D-5

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
     j0lmn=jlmn*(jlmn-1)/2
     do ilmn=1,jlmn
       il=indlmn(1,ilmn)
       im=indlmn(2,ilmn)
       iln=indlmn(5,ilmn)
       ilm=indlmn(4,ilmn)
       klmn=j0lmn+ilmn
       kln = pawtab%indklmn(2,klmn)

  !    Computation of (<phi_i|phi_j>-<tphi_i|tphi_j>)/r^3 radial integral

       ff(2:mesh_size)=(pawtab%phiphj(2:mesh_size,kln)-&
  &     pawtab%tphitphj(2:mesh_size,kln))/pawrad%rad(2:mesh_size)**3
       call pawrad_deducer0(ff,mesh_size,pawrad)
       call simp_gen(intgr3,ff,pawrad)

       do idir = 1, 3

! matrix element <S il im|L_idir|S jl jm>
         call slxyzs(il,im,idir,jl,jm,lms)

         dijnd(2*klmn-1,1) = dijnd(2*klmn-1,1) + intgr3*dreal(lms)*nucdipmom(idir)*permeability
         dijnd(2*klmn,1) = dijnd(2*klmn,1) + intgr3*dimag(lms)*nucdipmom(idir)*permeability

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
!!  cplex_dij=1 if dij is REAL, 2 if complex (2 for spin-orbit)
!!  ndij= number of spin components for Dij^SO
!!  nspden=number of spin density components
!!  paw_an <type(paw_an_type)>=paw arrays given on angular mesh, for current atom
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data, for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data, for current atom
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  spnorbscl=scaling factor for spin-orbit coupling
!!  vh1(mesh_size,v_size,nspden)=all-electron on-site Hartree potential for current atom
!!                     only spherical moment is used
!!  vxc1(mesh_size,v_size,nspden)=all-electron on-site XC potential for current atom
!!                                given on a (r,theta,phi) grid (v_size=angl_size)
!!                                or on (l,m) spherical moments (v_size=lm_size)
!!
!! OUTPUT
!!  dijso(cplex_dij*lmn2_size,ndij)= spin-orbit Dij terms
!!  cplex_dij=2 must be 2
!!        dijso(:,:,1) contains Dij_SO^up-up
!!        dijso(:,:,2) contains Dij_SO^dn-dn
!!        dijso(:,:,3) contains Dij_SO^up-dn
!!        dijso(:,:,4) contains Dij_SO^dn-up
!!
!! PARENTS
!!      m_pawdij,pawdenpot
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijso(cplex_dij,dijso,ndij,nspden,&
&                   pawang,pawrad,pawtab,pawxcdev,spnorbscl,vh1,vxc1)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdijso'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex_dij,ndij,nspden,pawxcdev
 real(dp), intent(in) :: spnorbscl
 type(pawang_type),intent(in) :: pawang
!arrays
 real(dp),intent(out) :: dijso(:,:)
 real(dp),intent(in) :: vh1(:,:,:),vxc1(:,:,:)
 type(pawrad_type),intent(in) :: pawrad
 type(pawtab_type),target,intent(in) :: pawtab
!Local variables ---------------------------------------
!scalars
 integer :: angl_size,cplex,idij,ij_size,ilm,ipts,ispden,jlm,klm,klmn,klmn1,kln
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
 cplex=1 ! DijSO exists only for GS
 nsploop=4

!Check data consistency
 if (cplex_dij/=2) then
   msg='cplex_dij must be 2 for spin-orbit coupling !'
   MSG_BUG(msg)
 end if
 if (ndij/=4) then
   msg='ndij must be 4 for spin-orbit coupling !'
   MSG_BUG(msg)
 end if
 if (pawang%use_ls_ylm==0) then
   msg='pawang%use_ls_ylm should be /=0 !'
   MSG_BUG(msg)
 end if
 if (size(dijso,1)/=cplex_dij*lmn2_size.or.size(dijso,2)/=ndij) then
   msg='invalid sizes for DijSO !'
   MSG_BUG(msg)
 end if
 if (size(vh1,1)/=cplex*mesh_size.or.size(vh1,2)<1.or.size(vh1,3)<1) then
   msg='invalid sizes for vh1 !'
   MSG_BUG(msg)
 end if
 if (size(vxc1,1)/=cplex*mesh_size.or.size(vxc1,3)/=nspden.or.&
&   (size(vxc1,2)/=angl_size.and.pawxcdev==0).or.&
&   (size(vxc1,2)/=lm_size.and.pawxcdev/=0)) then
   msg='invalid sizes for vxc1 !'
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
!!  cplex_dij=1 if dij is REAL, 2 if complex (2 for spin-orbit)
!!  ndij= number of spin components for Dij^hat
!!  nspden=number of spin density components
!!  nsppol=number of independent spin WF components
!!  pawtab <type(pawtab_type)>=paw tabulated starting data, for current atom
!!  vpawu(cplex_dij,lpawu*2+1,lpawu*2+1,nspden)=moments of LDA+U potential for current atom
!!  --- Optional arguments ---
!!    atvshift(natvshift,nsppol)=potential energy shift for lm channel & spin (current atom)
!!    fatvshift=factor that multiplies atvshift
!!    natvshift=number of atomic potential energy shifts (per atom)
!! 
!! OUTPUT
!!  dijpawu(cplex_dij*lmn2_size,ndij)=  D_ij^XC terms
!!   
!! PARENTS
!!      m_pawdij
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdiju(cplex_dij,dijpawu,ndij,nspden,nsppol,pawtab,vpawu,&
&                  natvshift,atvshift,fatvshift) ! optional arguments


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdiju'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex_dij,ndij,nspden,nsppol
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
 if (size(dijpawu,1)/=cplex_dij*lmn2_size.or.size(dijpawu,2)/=ndij) then
   msg='invalid sizes for dijpawu !'
   MSG_BUG(msg)
 end if
 if (size(vpawu,1)/=cplex_dij.or.size(vpawu,2)/=2*lpawu+1.or.&
&    size(vpawu,3)/=2*lpawu+1.or.size(vpawu,4)/=nspden) then
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
             dijpawu_idij(klmn1)=pawtab%phiphjint(icount)*coeffpawu(1) ! *dtset%userra
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

   if (nspden==4.or.ndij==4.or.cplex_dij==2) then
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
!!  cplex=(RF calculations only) - 1 if RF 1st-order quantities are REAL, 2 if COMPLEX
!!  cplex_dij=1 if dij is REAL, 2 if complex (2 for spin-orbit)
!!  lmselect(lm_size)=select the non-zero LM-moments of on-site potentials
!!  ndij= number of spin components for Dij^hat
!!  nspden=number of spin density components
!!  nsppol=number of independent spin WF components
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data, for current atom
!!  pawtab <type(pawtab_type)>=paw tabulated starting data, for current atom
!!  vpawx(1,lmn2_size,nspden)=moments of exact exchange potential
!!                    for current atom and for correlated electrons
!!  vxc_ex(mesh_size,lm_size,nspden)=all-electron on-site XC potential for current atom
!!                    taken into account only valence correlated electrons
!!
!! OUTPUT
!!  dijexxc(cplex_dij*lmn2_size,ndij)=  D_ij^XC terms
!!
!! NOTES
!!  cplex is for RF, cplex_dij is for non-collinear (nspinor==2)
!!
!! PARENTS
!!      m_pawdij
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijexxc(cplex,cplex_dij,dijexxc,lmselect,ndij,nspden,nsppol,&
&                      pawang,pawrad,pawtab,vpawx,vxc_ex)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdijexxc'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,cplex_dij,ndij,nspden,nsppol
 type(pawang_type),intent(in) :: pawang
!arrays
 logical :: lmselect(:)
 real(dp),intent(in) :: vpawx(:,:,:),vxc_ex(:,:,:)
 real(dp),intent(out) :: dijexxc(:,:)
 type(pawrad_type),intent(in) :: pawrad
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: icount,idij,idijend,ij_size,iln,in1,in2,ir,ir1,isel,ispden
 integer :: jln,j0ln,klm,klm1,klmn,klmn1,kln,lexexch,ln_min,ln_max,lmax,lmin
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
 if (size(dijexxc,1)/=cplex_dij*lmn2_size.or.size(dijexxc,2)/=ndij) then
   msg='invalid sizes for dijexxc !'
   MSG_BUG(msg)
 end if
 if (size(lmselect)/=lm_size) then
   msg='invalid size for lmselect !'
   MSG_BUG(msg)
 end if
 if (size(vxc_ex,1)/=cplex*mesh_size.or.size(vxc_ex,2)/=lm_size.or.&
&    size(vxc_ex,3)/=nspden) then
   msg='invalid sizes for vxc_ex !'
   MSG_BUG(msg)
 end if
 if (size(vpawx,1)/=1.or.size(vpawx,2)/=lmn2_size.or.&
&    size(vpawx,3)/=nspden) then
   msg='invalid sizes for vpawx !'
   MSG_BUG(msg)
 end if
 if (cplex_dij<cplex) then
   msg='cplex_dij must be >= cplex !'
   MSG_BUG(msg)
 end if

!Init memory
 dijexxc=zero
 LIBPAW_ALLOCATE(dijexxc_idij,(cplex*lmn2_size))
 LIBPAW_ALLOCATE(vxcij1,(cplex*ij_size))
 LIBPAW_ALLOCATE(ff,(mesh_size))
 LIBPAW_ALLOCATE(gg,(mesh_size))

!----------------------------------------------------------
!Loop over spin components
!----------------------------------------------------------
 nsploop=nsppol;if (ndij==4) nsploop=4
 do idij=1,nsploop

   if (idij<=nsppol.or.(nspden==4.and.idij<=3).or.cplex==2) then

     idijend=idij+idij/3;if (cplex==2) idijend=idij
     do ispden=idij,idijend

       dijexxc_idij=zero

!      ----------------------------------------------------------
!      Summing over (l,m) moments
!      ----------------------------------------------------------
       do klm=1,lm_size
         if (lmselect(klm)) then

!          ===== Vxc_ij_1 (tmp) =====
           vxcij1=zero
           if (cplex==1) then
             do jln=ln_min,ln_max
               j0ln=jln*(jln-1)/2
               do iln=ln_min,jln
                 kln=j0ln+iln
                 ff(1:mesh_size)= &
&                  vxc_ex(1:mesh_size,klm,idij)*pawtab%phiphj(1:mesh_size,kln)
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
&                   vxc_ex(ir1-1,klm,ispden)*pawtab%phiphj(ir,kln)
                   gg(ir)= &
&                   vxc_ex(ir1,klm,ispden)*pawtab%phiphj(ir,kln)
                 end do
                 call simp_gen(vxcij1(2*kln-1),ff,pawrad)
                 call simp_gen(vxcij1(2*kln  ),gg,pawrad)
               end do
             end do
           end if

!          ===== Accumulate Vxc_ij_1 over klm moments =====
           if (cplex==1) then
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
           else ! cplex==2
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
               klmn1=klmn1+cplex
             end do ! Loop klmn
           end if

         end if ! lmselect
       end do  ! Loop klm

!      Mix Hartree and GGA terms
       if (cplex==1) then
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
       else ! cplex=2
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
           klmn1=klmn1+cplex
         end do ! Loop klmn
       end if

!      ----------------------------------------------------------
!      Deduce some part of Dij according to symmetries
!      ----------------------------------------------------------

       if (cplex==1) then
         if (ispden<3) then
           if (cplex_dij==1) then
             dijexxc(1:lmn2_size,idij)=dijexxc_idij(1:lmn2_size)
           else
             klmn1=1
             do klmn=1,lmn2_size
               dijexxc(klmn1  ,idij)=dijexxc_idij(klmn)
               dijexxc(klmn1+1,idij)=zero
               klmn1=klmn1+cplex_dij
             end do
           end if
         else
           klmn1=max(1,ispden-2)
           do klmn=1,lmn2_size
             dijexxc(klmn1,idij)=dijexxc_idij(klmn)
             klmn1=klmn1+cplex_dij
           end do
         end if
       else !cplex=2
         if (ispden<=3) then
           dijexxc(1:cplex*lmn2_size,idij)=dijexxc_idij(1:cplex*lmn2_size)
         else     
           klmn1=1  ! Remember V(4) contains i.V^21
           do klmn=1,lmn2_size
             dijexxc(klmn1  ,idij)= dijexxc_idij(klmn+1)
             dijexxc(klmn1+1,idij)=-dijexxc_idij(klmn  )
             klmn1=klmn1+cplex_dij
           end do
         end if
       end if

     end do !ispden

   else if (nspden==4.and.idij==4) then ! cplex=1 here
     dijexxc(:,idij)=dijexxc(:,idij-1)
     klmn1=2
     do klmn=1,lmn2_size
       dijexxc(klmn1,idij)=-dijexxc(klmn1,idij)
       klmn1=klmn1+cplex_dij
     end do
   else if (nsppol==1.and.idij==2) then ! cplex=1 here
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
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL; if 2, COMPLEX
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!  idir=direction of atomic displacement (in case of phonons perturb.)
!!  ipert=nindex of perturbation
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  mpi_comm_grid=--optional-- MPI communicator over real space grid components
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms
!!  option=0: computes full frozen part of Dij
!!         1: computes frozen part of Dij without contribution from Vpsp1
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(my_natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  qphon(3)=wavevector of the phonon
!!  rprimd(3,3)=dimensional primitive translations for real space
!!  ucvol=unit cell volume (bohr^3)
!!  vpsp1(cplex*nfft)= first-order change of local potential
!!  vtrial(nfft,nspden)= total GS potential
!!  vxc(nfft,nspden)=XC potential
!!  xred(3,my_natom)= reduced atomic coordinates
!!
!! OUTPUT
!!  paw_ij1(iatom)%dijfr(cplex_dij*lmn2_size,nspden)=
!!                  frozen contribution to psp strength Dij
!!                  =Int_R^3{vtrial*Sum_LM[Q_ij_q^LM^(1)] + Vloc^(1)*Sum_LM[Q_ij_q^LM]}
!!
!! PARENTS
!!      d2frnl,dfpt_nstpaw,dfpt_rhofermi,dfpt_scfcv
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

subroutine pawdijfr(cplex,gprimd,idir,ipert,my_natom,natom,nfft,ngfft,nspden,ntypat,&
&          option,paw_ij1,pawang,pawfgrtab,pawrad,pawtab,qphon,rprimd,ucvol,vpsp1,vtrial,vxc,xred,&
&          mpi_atmtab,comm_atom,mpi_comm_grid) ! optional arguments (parallelism)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdijfr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,my_natom,natom,nfft,nspden,ntypat,option
 integer,optional,intent(in) :: comm_atom,mpi_comm_grid
 real(dp),intent(in) :: ucvol
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gprimd(3,3),qphon(3),rprimd(3,3)
 real(dp),intent(in) :: vpsp1(cplex*nfft),vtrial(nfft,nspden),vxc(nfft,nspden)
 real(dp),intent(in) :: xred(3,natom)
 type(paw_ij_type),intent(inout) :: paw_ij1(my_natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: cplex_dij,dplex,iatom,iatom_tot,ic,ier,ils,ilslm,isel,ispden,istr,itypat,jc
 integer :: klm,klmn,klmn1,kln,lm_size,lmn2_size,lm0,lmax,lmin,mesh_size
 integer :: mm,my_comm_atom,my_comm_grid,mu,mua,mub,nfftot,nfgd,optgr0,optgr1,optgr2,usexcnhat
 logical :: has_phase,my_atmtab_allocated,need_dijfr_1,need_dijfr_2,need_dijfr_3,need_dijfr_4
 logical :: paral_atom,qne0,testdij1,testdij2,testdij3
 real(dp) :: c1,fact,intg,rg1
 character(len=500) :: msg
!arrays
 integer,parameter :: m_index(3)=(/1,-1,0/)
 integer,pointer :: my_atmtab(:)
 integer,parameter :: alpha(9)=(/1,2,3,3,3,2,2,1,1/),beta(9)=(/1,2,3,2,1,1,3,3,2/)
 real(dp) :: contrib(2)
 real(dp),allocatable :: ff(:),intv(:,:),intvloc(:,:),intv_tmp(:,:),rg(:),vloc(:,:)

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
   if (paw_ij1(1)%cplex_dij<cplex) then
     msg='paw_ij1()%cplex_dij must be >=cplex !'
     MSG_BUG(msg)
   end if
   if (paw_ij1(1)%cplex/=cplex) then
     msg='paw_ij1()%cplex and cplex must be equal !'
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
 dplex=cplex-1

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
   has_phase=.false.
   if (need_dijfr_2) then
     if (qne0.and.(pawfgrtab(iatom)%expiqr_allocated==0)) then
       if (allocated(pawfgrtab(iatom)%expiqr))  then
         LIBPAW_DEALLOCATE(pawfgrtab(iatom)%expiqr)
       end if
       LIBPAW_ALLOCATE(pawfgrtab(iatom)%expiqr,(2,nfgd))
       call pawexpiqr(pawfgrtab(iatom)%expiqr,gprimd,nfgd,qphon,&
&                     pawfgrtab(iatom)%rfgd,xred(:,iatom_tot))
       pawfgrtab(iatom)%expiqr_allocated=2
     end if
     has_phase=(pawfgrtab(iatom)%expiqr_allocated/=0)
   end if

!  Loop over spin components
   do ispden=1,nspden

!    ============ Phonons ====================================
     if (ipert<=natom) then

       if (need_dijfr_1.or.need_dijfr_2) then

         LIBPAW_ALLOCATE(intv,(cplex,lm_size))
         intv(:,:)=zero

!        First part: Int_R^3{vtrial*Sum_LM[Q_ij_q^LM^(1)]}
         if (need_dijfr_1) then

!          ----- Retrieve potential Vlocal (subtle if nspden=4 ;-)
           if (nspden/=4) then
             LIBPAW_ALLOCATE(vloc,(1,nfgd))
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
             LIBPAW_ALLOCATE(vloc,(2,nfgd))
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

!          ----- Compute Integral [ Vtrial(r).(g_l(r).Y_lm(r))^(1) dr ]
           LIBPAW_ALLOCATE(intv_tmp,(cplex,3))
           do ilslm=1,lm_size
             intv_tmp=zero
             if (nspden/=4) then
               do ic=1,nfgd
                 do mu=1,3
!                  Minus sign because dg(r-R)/dR = -dg(r-R)/dr
                   contrib(1)=-vloc(1,ic)*pawfgrtab(iatom)%gylmgr(mu,ic,ilslm)
                   intv_tmp(1,mu)=intv_tmp(1,mu)+contrib(1)
                 end do
               end do
             else ! nspden=4
               do ic=1,nfgd
                 do mu=1,3
!                  Minus sign because dg(r-R)/dR = -dg(r-R)/dr
                   contrib(1:2)=-vloc(1:2,ic)*pawfgrtab(iatom)%gylmgr(mu,ic,ilslm)
                   intv_tmp(1:2,mu)=intv_tmp(1:2,mu)+contrib(1:2)
                 end do
               end do
             end if
!            Convert from cartesian to reduced coordinates
             intv(1:cplex,ilslm)=intv(1:cplex,ilslm) &
&             +(rprimd(1,idir)*intv_tmp(1:cplex,1) &
&             +rprimd(2,idir)*intv_tmp(1:cplex,2) &
&             +rprimd(3,idir)*intv_tmp(1:cplex,3))
           end do
           LIBPAW_DEALLOCATE(vloc)
           LIBPAW_DEALLOCATE(intv_tmp)
         end if ! need_dijfr_1

!        2nd part: Int_R^3{Vloc^(1)*Sum_LM[Q_ij_q^LM]}
         if (need_dijfr_2) then

           if (ispden==1) then

!            ----- Retrieve potential Vloc^(1)
             LIBPAW_ALLOCATE(vloc,(cplex,nfgd))
             do ic=1,nfgd
               jc=cplex*pawfgrtab(iatom)%ifftsph(ic)-dplex
               vloc(1:cplex,ic)=vpsp1(jc:jc+dplex)
             end do

!            ----- Compute Integral [ Vloc^(1)(r).g_l(r).Y_lm(r) ]
             LIBPAW_ALLOCATE(intvloc,(cplex,lm_size))
             intvloc=zero
             if (has_phase) then
               if (cplex==1) then
                 do ilslm=1,lm_size
                   do ic=1,nfgd
                     contrib(1)=vloc(1,ic)*pawfgrtab(iatom)%gylm(ic,ilslm)
                     intvloc(1,ilslm)=intvloc(1,ilslm)+contrib(1)*pawfgrtab(iatom)%expiqr(1,ic)
                   end do
                 end do
               else
                 do ilslm=1,lm_size
                   do ic=1,nfgd
                     contrib(1:cplex)=vloc(1:cplex,ic)*pawfgrtab(iatom)%gylm(ic,ilslm)
                     intvloc(1,ilslm)=intvloc(1,ilslm)+contrib(1)*pawfgrtab(iatom)%expiqr(1,ic) &
&                     -contrib(2)*pawfgrtab(iatom)%expiqr(2,ic)
                     intvloc(2,ilslm)=intvloc(2,ilslm)+contrib(1)*pawfgrtab(iatom)%expiqr(2,ic) &
&                     +contrib(2)*pawfgrtab(iatom)%expiqr(1,ic)
                   end do
                 end do
               end if
             else ! no phase
               do ilslm=1,lm_size
                 do ic=1,nfgd
                   contrib(1:cplex)=vloc(1:cplex,ic)*pawfgrtab(iatom)%gylm(ic,ilslm)
                   intvloc(1:cplex,ilslm)=intvloc(1:cplex,ilslm)+contrib(1:cplex)
                 end do
               end do
             end if
             LIBPAW_DEALLOCATE(vloc)
           end if ! ispden=1

           if (ispden<=min(nspden,2)) then
             intv(1:cplex,1:lm_size)=intv(1:cplex,1:lm_size)+intvloc(1:cplex,1:lm_size)
             if (ispden==min(nspden,2))  then
               LIBPAW_DEALLOCATE(intvloc)
             end if
           end if
         end if ! need_dijfr_2

!        Apply ucvol/nfft factor on integral
         intv(:,:)=fact*intv(:,:)

!        --- Reduction in case of parallelization ---
         call xmpi_sum(intv,my_comm_grid,ier)

         paw_ij1(iatom)%dijfr(:,ispden)=zero

!        ---- Loop over (i,j) components
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
!                The following works only because cplex<=cplex_dij
!                if cplex<cplex_dij, add zero to imaginary par of dijfr
                 paw_ij1(iatom)%dijfr(klmn1:klmn1+dplex,ispden)= &
&                 paw_ij1(iatom)%dijfr(klmn1:klmn1+dplex,ispden) &
&                 +pawtab(itypat)%qijl(ilslm,klmn)*intv(1:cplex,ilslm)
               end if
             end do
           end do
           klmn1=klmn1+cplex_dij
         end do
         LIBPAW_DEALLOCATE(intv)

!        Dijfr is marked as computed
         paw_ij1(iatom)%has_dijfr=2
       end if

!    ============ Electric field perturbation =======================
     else if (ipert==natom+2.or.ipert==natom+11) then

       if (need_dijfr_3) then

!        The following factor arises in expanding the angular dependence of the dipole
!        vector in terms of real spherical harmonics. The real spherical harmonics are as
!        in the routine initylmr.F90;
!        see http://www.unioviedo.es/qcg/art/Theochem419-19-ov-BF97-rotation-matrices.pdf
         c1 = sqrt(four_pi/three)
         mesh_size=pawtab(itypat)%mesh_size

         if (ispden==1) then

           LIBPAW_ALLOCATE(ff,(mesh_size))
           LIBPAW_ALLOCATE(rg,(3))

!          loop over basis state pairs for this atom
           klmn1=1
           do klmn = 1, paw_ij1(iatom)%lmn2_size
             klm =pawtab(itypat)%indklmn(1,klmn)
             kln =pawtab(itypat)%indklmn(2,klmn)
             lmin=pawtab(itypat)%indklmn(3,klmn)
             lmax=pawtab(itypat)%indklmn(4,klmn)

!            Select only l=1, because the dipole is a vector operator
             if (lmin==1) then
               lm0=3  ! (l^2+l+1) for l=1

!              Computation of <phi_i|r|phi_j>- <tphi_i|r|tphi_j>
!              the dipole vector has radial dependence r
               ff(1:mesh_size)=(pawtab(itypat)%phiphj(1:mesh_size,kln)&
&               -pawtab(itypat)%tphitphj(1:mesh_size,kln))&
&               *pawrad(itypat)%rad(1:mesh_size)
!               call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
               call simp_gen(intg,ff,pawrad(itypat))

!              Compute <S_li_mi|r-R|S_lj_mj>: use a real Gaunt expression (with selection rule)
               rg(1:3)=zero
               do ic=1,3
                 isel=pawang%gntselect(lm0+m_index(ic),klm)
                 if (isel>0) rg(ic)=pawang%realgnt(isel)
               end do

!              Translate from cartesian to reduced coordinates (in idir direction)
               rg1=gprimd(1,idir)*rg(1)+gprimd(2,idir)*rg(2)+gprimd(3,idir)*rg(3)

!              Build sqrt(4pi/3).<S_li_mi|r-R|S_lj_mj>.(<phi_i|r-R|phi_j>- <tphi_i|r-R|tphi_j>
               paw_ij1(iatom)%dijfr(klmn1,ispden)=c1*rg1*intg
               if (cplex_dij==2) paw_ij1(iatom)%dijfr(klmn1+1,ispden)=zero

             else
               paw_ij1(iatom)%dijfr(klmn1,ispden)=zero
             end if ! end gaunt constraint

             klmn1=klmn1+cplex_dij
           end do ! end loop over lmn2_size pairs of basis states
           LIBPAW_DEALLOCATE(ff)
           LIBPAW_DEALLOCATE(rg)

!          Dijfr is spin-independent for electric field case
         else if (ispden==2) then
           paw_ij1(iatom)%dijfr(:,ispden)=paw_ij1(iatom)%dijfr(:,1)
         else
           paw_ij1(iatom)%dijfr(:,ispden)=zero
         end if

!        Dijfr is marked as computed
         paw_ij1(iatom)%has_dijfr=2
       end if

!    ============ Elastic tensor ===============================
     else if (ipert==natom+3.or.ipert==natom+4) then
!     ----- Retrieve potential Vlocal (subtle if nspden=4 ;-)
       if (nspden/=4) then
         LIBPAW_ALLOCATE(vloc,(1,nfgd))
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
         LIBPAW_ALLOCATE(vloc,(2,nfgd))
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

       LIBPAW_ALLOCATE(intv,(cplex,lm_size))
       intv(:,:)=zero
!      option = 0 Insulator case
       if(option==0)then
         do ilslm=1,lm_size
           do ic=1,nfgd
             jc=pawfgrtab(iatom)%ifftsph(ic) 
             contrib(1) = 0
!            Int_R^3{vtrial*Sum_LM[Q_ij_q^LM^(1)]}
             mua=alpha(istr);mub=beta(istr)
             contrib(1)=contrib(1)+half*vloc(1,ic)&
&              *(pawfgrtab(iatom)%gylmgr(mua,ic,ilslm)*pawfgrtab(iatom)%rfgd(mub,ic)&
&              + pawfgrtab(iatom)%gylmgr(mub,ic,ilslm)*pawfgrtab(iatom)%rfgd(mua,ic))

!            Int_R^3{Vloc^(1)*Sum_LM[Q_ij_q^LM]}
             contrib(1)=contrib(1)+vpsp1(jc)*pawfgrtab(iatom)%gylm(ic,ilslm)

!            delta_{alphabeta}Int_R^3{Vloc*Sum_LM[Q_ij_q^LM]}
             if(istr<=3)then
               contrib(1)=contrib(1)+vloc(1,ic)*pawfgrtab(iatom)%gylm(ic,ilslm)
             end if

             intv(1,ilslm)=intv(1,ilslm)+contrib(1)
           end do
         end do

!      option = 1 Metal case (without Vpsp1)
       else if (option==1)then
         do ilslm=1,lm_size
           do ic=1,nfgd
             jc=pawfgrtab(iatom)%ifftsph(ic) 
             contrib(1) = 0
!            Int_R^3{vtrial*Sum_LM[Q_ij_q^LM^(1)]}
             mua=alpha(istr);mub=beta(istr)
             contrib(1)=contrib(1)+half*vloc(1,ic)&
&              *(pawfgrtab(iatom)%gylmgr(mua,ic,ilslm)*pawfgrtab(iatom)%rfgd(mub,ic)&
&              + pawfgrtab(iatom)%gylmgr(mub,ic,ilslm)*pawfgrtab(iatom)%rfgd(mua,ic))

!            delta_{alphabeta}Int_R^3{Vtrial*Sum_LM[Q_ij_q^LM]}
             if(istr<=3)then
               contrib(1)=contrib(1)+vloc(1,ic)*pawfgrtab(iatom)%gylm(ic,ilslm)
             end if

             intv(1,ilslm)=intv(1,ilslm)+contrib(1)
           end do
         end do
       end if

!      Apply ucvol/nfft factor on integral
       intv(:,:)=fact*intv(:,:)

!      --- Reduction in case of parallelization ---
       call xmpi_sum(intv,my_comm_grid,ier)

       paw_ij1(iatom)%dijfr(:,ispden)=zero

!      ---- Loop over (i,j) components
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
!              The following works only because cplex<=cplex_dij
!              if cplex<cplex_dij, add zero to imaginary par of dijfr
               paw_ij1(iatom)%dijfr(klmn1:klmn1+dplex,ispden)= &
&                paw_ij1(iatom)%dijfr(klmn1:klmn1+dplex,ispden) &
&                +pawtab(itypat)%qijl(ilslm,klmn)*intv(1:cplex,ilslm)
             end if
           end do
         end do
         klmn1=klmn1+cplex_dij
       end do
       LIBPAW_DEALLOCATE(intv)
       LIBPAW_DEALLOCATE(vloc)
!      Dijfr is marked as computed
       paw_ij1(iatom)%has_dijfr=2

     end if ! ipert

!    End loop over spin components
   end do ! ispden

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
!!  cplex_dij=1 if dij is REAL, 2 if complex (2 for spin-orbit)
!!  ndij=number of spin components for Dij^hat
!!  nspden=number of spin density components
!!  pawprtvol=control print volume and debugging output for PAW
!!  noccmmp(cplex_dij,2*lpawu+1,2*lpawu+1,ndij)=density matrix in the augm. region
!!  nocctot(ndij)=number of electrons in the correlated subspace
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data:
!!     %usepawu, %upawu, %jpau
!!     %vee(2*lpawu+1*4)=screened coulomb matrix
!!
!! OUTPUT
!!  vpawu(cplex_dij,lpawu*2+1,lpawu*2+1,nspden)=lda+u potential
!!                                 (see eg PRB 52, 5467 (1995))
!!
!! PARENTS
!!      ldau_self,m_pawdij,m_pawhr
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_allgatherv
!!
!! SOURCE

 subroutine pawpupot(cplex_dij,ndij,noccmmp,nocctot,nspden,&
&                    pawprtvol,pawtab,vpawu)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawpupot'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex_dij,ndij,nspden,pawprtvol
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

 integer :: iplex,ispden,jspden,lpawu,m1,m11,m2,m21,m3,m31,m4,m41
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
 if(option_interaction==3.and.pawtab%usepawu>=10) then
   msg = "Option_interaction==3 is not compatible with usepawu>=10 in pawpupot"
   MSG_ERROR(msg)
 end if
 if (size(vpawu,1)/=cplex_dij.or.size(vpawu,2)/=2*lpawu+1.or.&
&    size(vpawu,3)/=2*lpawu+1.or.size(vpawu,4)/=nspden) then
   write (msg,'(a,4I5, a,4I5)') ' invalid sizes for vpawu !',cplex_dij,2*lpawu+1,2*lpawu+1,nspden, &
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
!cf PRB 52 5467 (1995)
!-----------------------------------------------------

 vpawu=zero
 do ispden=1,nspden

   if (ispden<=2) then   ! cases ndij=4, ispden=1,2 or ndij<4
     jspden=min(ndij,2)-ispden+1   ! (ispden,ndij)=(1,4)=>jspden=2

     if (nspden<=2) then
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
       if(pawtab%usepawu==1) then ! not activated if usepawu=10 !!
!        Here we compute vpawu=vpawu-v_dc 
         vpawu(1,m11,m11,ispden)=vpawu(1,m11,m11,ispden)-pawtab%upawu*(n_tot-half)
         if (ndij/=4.or.option_interaction==2) then
           vpawu(1,m11,m11,ispden)=vpawu(1,m11,m11,ispden)+pawtab%jpawu*(n_sig-half)
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

!  Non-collinear magnetism: add non-diagonal term; see (Eq 6) in PRB 72, 024458 (2005)
!  BA Here, we compute the transpose --- with respect to spin indices --- of
!  BA equation (6) of this reference, because of differences in notations,
!  BA namely Eband=\sum rhoij^{alpha,beta}*Dij(beta,alpha) contrary to PRB 72, 024458 (2005)
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
       if(pawtab%usepawu==1.and.option_interaction==3) then ! not activated if usepawu=10 !!
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
     if (ispden==nspden) then
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
!!  nspden=number of spin-density components
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

 subroutine pawxpot(nspden,pawprtvol,pawrhoij,pawtab,vpawx)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawxpot'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: nspden,pawprtvol
 type(pawrhoij_type),intent(in) :: pawrhoij
 type(pawtab_type),intent(in) :: pawtab
 real(dp),intent(out) :: vpawx(:,:,:)

!Local variables ---------------------------------------
!scalars
 integer :: irhoij,irhoij1,ispden,jrhoij,jrhoij1,klmn,klmn1,lexexch,ll,lmn2_size
 integer :: m11,m21,m31,m41,n1,n2,n3,n4,nk,nn1,nn2
 real(dp) :: tot
 character(len=500) :: msg
!arrays
 integer :: indn(3,3)
 real(dp) :: factnk(6)

! *****************************************************

!Useful data
 lexexch=pawtab%lexexch
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
&    size(vpawx,3)/=nspden) then
   msg='invalid sizes for vpawx !'
   MSG_BUG(msg)
 end if

!=====================================================
!Compute local exact exchange Potential
!on the basis of projectors.
!-----------------------------------------------------

 vpawx=zero
 do ispden=1,nspden
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
         jrhoij1=jrhoij1+pawrhoij%cplex
       end do !irhoij1
     end if
     jrhoij=jrhoij+pawrhoij%cplex
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
       jrhoij=jrhoij+pawrhoij%cplex
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
!!  paw_ij(natom)%cplex_dij=1 if dij are REAL, 2 if they are COMPLEX
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symdij'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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
 integer :: at_indx,cplex,cplex_dij,iafm,iatom,iatom_tot
 integer :: il,il0,ilmn,iln,iln0,ilpm,indexi,indexii,indexj,indexjj,indexjj0,indexk,indexkc
 integer :: iplex,irot,ispden,itypat,j0lmn,jl,jl0,jlmn,jln,jln0,jlpm,jspden
 integer :: klmn,klmnc,kspden,lmn_size,mi,mj,my_comm_atom,mu,natinc,ndij0,ndij1,nu,optsym,sz1,sz2
 logical,parameter :: afm_noncoll=.true.  ! TRUE if antiferro symmetries are used with non-collinear magnetism
 logical :: antiferro,have_phase,my_atmtab_allocated,noncoll,paral_atom,use_afm
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
 integer :: idum(0)
 real(dp) :: dijc(2),factsym(2),phase(2),rotdij(2,2),sumdij(2,2)
 real(dp),allocatable :: dijnew(:,:),dijtmp(:,:),rotmag(:,:),summag(:,:)
 real(dp),allocatable :: symrec_cart(:,:,:),work(:,:)
 character(len=7),parameter :: dspin(6)=(/"up     ","down   ","up-up  ", &
&                                         "dwn-dwn","up-dwn ","dwn-up "/)
 type(coeff2_type),target, allocatable :: my_tmp_dij(:)
 type(coeff2_type),pointer :: tmp_dij(:)

!DEBUG_ALTERNATE_ALGO
!integer :: i1,i2,i3,i4,symrel_conv(3,3)
!real(dp) :: spinrot(4)
!real(dp),allocatable :: dijtemp(:,:),sumrhoso(:,:)
!complex(dpc) :: dijt(2,2),dijt2(2,2),Rspinrot(2,2)
!DEBUG_ALTERNATE_ALGO

! *********************************************************************

!Tests of compatibility:
 if (my_natom>0) then
   if ((option_dij==1.and.paw_ij(1)%has_dijhat==0).or.&
&   (option_dij==2.and.paw_ij(1)%has_dijU==0).or.&
&   (option_dij==3.and.paw_ij(1)%has_dijxc==0).or.&
&   (option_dij==4.and.paw_ij(1)%has_dijxc_hat==0).or.&
&   (option_dij==5.and.paw_ij(1)%has_dijxc_val==0).or.&
&   (option_dij==6.and.paw_ij(1)%has_dijso==0).or.&
&   (option_dij==7.and.paw_ij(1)%has_dijexxc==0).or.&
&   (option_dij==8.and.paw_ij(1)%has_dijfr==0).or.&
&   (option_dij==9.and.paw_ij(1)%has_dijnd==0)) then
     msg='Incompatibilty between option_dij and allocation of Dij!'
     MSG_BUG(msg)
   end if
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Symmetrization occurs only when nsym>1
 if (nsym>1.and.ipert/=natom+1.and.ipert/=natom+10) then

   if (pawang%nsym==0) then
     msg='pawang%zarot must be allocated!'
     MSG_BUG(msg)
   end if

   cplex_dij=1;antiferro=.false.;noncoll=.false.
   if (my_natom>0) then
     cplex_dij=paw_ij(1)%cplex_dij
!    Antiferro case ?
     antiferro=(paw_ij(1)%nspden==2.and.paw_ij(1)%nsppol==1.and.paw_ij(1)%ndij/=4)
!    Non-collinear case
     noncoll=(paw_ij(1)%ndij==4)
     if (noncoll.and.paw_ij(1)%cplex_dij/=2) then
       msg='cplex_dij must be 2 with ndij=4!'
       MSG_BUG(msg)
     end if
   end if
   if (noncoll) then
     LIBPAW_ALLOCATE(summag,(cplex_dij,3))
     LIBPAW_ALLOCATE(rotmag,(cplex_dij,3))
     LIBPAW_ALLOCATE(work,(cplex_dij,3))
     LIBPAW_ALLOCATE(symrec_cart,(3,3,nsym))
     do irot=1,nsym
       symrec_cart(:,:,irot)=symdij_symcart(gprimd,rprimd,symrec(:,:,irot))
     end do
!DEBUG_ALTERNATE_ALGO
!    if(lsymnew) then
!      LIBPAW_ALLOCATE(sumrhoso,(cplex_dij,4))
!    end if
!DEBUG_ALTERNATE_ALGO
   end if
!  Do we use antiferro symmetries ?
   use_afm=((antiferro).or.(noncoll.and.afm_noncoll))

!  Do we have a phase due to q-vector (phonons only) ?
   have_phase=.false.
   if (ipert>0.and.present(qphon).and.my_natom>0) then
     have_phase=(abs(qphon(1))>tol8.or.abs(qphon(2))>tol8.or.abs(qphon(3))>tol8)
     if (have_phase.and.cplex_dij==1) then
       msg='Should have cplex_dij=2 for a non-zero q!'
       MSG_BUG(msg)
     end if
   end if

!  Have to make a temporary copy of dij
   LIBPAW_DATATYPE_ALLOCATE(my_tmp_dij,(my_natom))
   do iatom=1,my_natom
     sz1=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size;sz2=paw_ij(iatom)%ndij
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
     end if
     !Has to translate Dij^{alpha,beta} into (Dij, Dij magnetic field) format
     if (paw_ij(1)%ndij==4) then
       my_tmp_dij(iatom)%value(:,1)=dijtmp(:,1)+dijtmp(:,2)
       my_tmp_dij(iatom)%value(:,2)=dijtmp(:,3)+dijtmp(:,4)
       my_tmp_dij(iatom)%value(:,4)=dijtmp(:,1)-dijtmp(:,2)
       do klmn=1,paw_ij(iatom)%lmn2_size
         my_tmp_dij(iatom)%value(2*klmn-1,3)=-dijtmp(2*klmn  ,3)+dijtmp(2*klmn  ,4)
         my_tmp_dij(iatom)%value(2*klmn  ,3)= dijtmp(2*klmn-1,3)-dijtmp(2*klmn-1,4)
       end do
!DEBUG_ALTERNATE_ALGO
!      if(lsymnew) my_tmp_dij(iatom)%value(:,:)=dijtmp(:,:)
!DEBUG_ALTERNATE_ALGO
     else
       my_tmp_dij(iatom)%value(:,:)=dijtmp(:,:)
     end if
     LIBPAW_DEALLOCATE(dijtmp)
   end do

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

   ndij1=1
   if (antiferro) ndij1=2
   if (noncoll)   ndij1=4
   ndij0=ndij1-1
   LIBPAW_ALLOCATE(dijnew,(cplex_dij,ndij1))

!  Loops over atoms and spin components
   do iatom=1,my_natom
     iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
     itypat=paw_ij(iatom)%itypat
     lmn_size=paw_ij(iatom)%lmn_size
     cplex_dij=paw_ij(iatom)%cplex_dij
     cplex=paw_ij(iatom)%cplex
     indlmn => pawtab(itypat)%indlmn

!DEBUG_ALTERNATE_ALGO
!    if (noncoll.and.lsymnew) then
!      LIBPAW_ALLOCATE(dijtemp,(paw_ij(iatom)%cplex_dij,paw_ij(iatom)%ndij))
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

           rotdij(:,:)=zero
           if (noncoll) rotmag(:,:)=zero
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

             if (have_phase) then
               arg=two_pi*(qphon(1)*indsym(1,irot,iatom)+qphon(2)*indsym(2,irot,iatom) &
&                         +qphon(3)*indsym(3,irot,iatom))
               phase(1)=cos(arg);phase(2)=sin(arg)
             end if

             sumdij(:,:)=zero
             if (noncoll) summag(:,:)=zero

!            Accumulate values over (mi,mj) and symmetries
             do mj=1,2*jl+1
               indexjj=indexj+mj;indexjj0=indexjj*(indexjj-1)/2
               do mi=1,2*il+1
                 indexii=indexi+mi
                 factsym(:)=one
                 if (indexii<=indexjj) then
                   indexk=indexjj0+indexii
                   if(cplex_dij==2.and.cplex==1) factsym(2)=one
                 else
                   indexk=indexii*(indexii-1)/2+indexjj
                   if(cplex_dij==2.and.cplex==1) factsym(2)=-one
                 end if
                 indexkc=cplex_dij*(indexk-1)
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
                   sumdij(1,iafm)=sumdij(1,iafm)+factsym(1)*zarot2*tmp_dij(at_indx)%value(indexkc+1,kspden)
                   if(cplex_dij==2) sumdij(2,iafm)= &
&                         sumdij(2,iafm)+factsym(2)*factafm*zarot2*tmp_dij(at_indx)%value(indexkc+2,kspden)
                 end if

                 if (noncoll.and.(.not.lsymnew)) then
                   do mu=1,3
                     summag(1,mu)=summag(1,mu)+factsym(1)*factafm*zarot2*tmp_dij(at_indx)%value(indexkc+1,1+mu)
                     if(cplex_dij==2) summag(2,mu)= &
&                           summag(2,mu)+factsym(2)*zarot2*tmp_dij(at_indx)%value(indexkc+2,1+mu)
                   end do
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
             if (have_phase) then
               if((.not.noncoll).or.(.not.lsymnew)) then
                 dijc(1:2)=sumdij(1:2,iafm)
                 sumdij(1,iafm)=phase(1)*dijc(1)-phase(2)*dijc(2)
                 sumdij(2,iafm)=phase(1)*dijc(2)+phase(2)*dijc(1)
               end if
               if (noncoll.and.(.not.lsymnew)) then
                 do mu=1,3
                   dijc(1:2)=summag(1:2,mu)
                   summag(1,mu)=phase(1)*dijc(1)-phase(2)*dijc(2)
                   summag(2,mu)=phase(1)*dijc(2)+phase(2)*dijc(1)
                 end do
               end if
!DEBUG_ALTERNATE_ALGO
!              if (noncoll.and.(lsymnew) then
!                do mu=1,4
!                  dijc(1:2)=sumrhoso(1:2,mu)
!                  sumrhoso(1,mu)=phase(1)*dijc(1)-phase(2)*dijc(2)
!                  sumrhoso(2,mu)=phase(1)*dijc(2)+phase(2)*dijc(1)
!                  end do
!                end do
!              end if
!DEBUG_ALTERNATE_ALGO
             end if

!            Add contribution of this rotation
             rotdij(1:cplex_dij,iafm)=rotdij(1:cplex_dij,iafm)+sumdij(1:cplex_dij,iafm)
             if (noncoll.and.(.not.lsymnew)) then
!              If non-collinear case, rotate Dij magnetization
!              Should use symrel^1 but use transpose[symrec] instead
               do nu=1,3
                 do mu=1,3
                   !we need the transpose ?
                   rotmag(1:cplex_dij,mu)=rotmag(1:cplex_dij,mu) &
&                       +symrec_cart(mu,nu,irot)*summag(1:cplex_dij,nu)
                 end do
               end do
             end if

           end do ! End loop over symmetries

           if((.not.noncoll).or.(.not.lsymnew)) then
!            Store new value of dij
             do iplex=1,cplex_dij
               dijnew(iplex,1)=rotdij(iplex,1)/nsym_used(1)
               if (abs(dijnew(iplex,1))<=tol10) dijnew(iplex,1)=zero
             end do

!            Antiferromagnetic case: has to fill up "down" component of dij
             if (antiferro.and.nsym_used(2)>0) then
               do iplex=1,cplex_dij
                 dijnew(iplex,2)=rotdij(iplex,2)/nsym_used(2)
                 if (abs(dijnew(iplex,2))<=tol10) dijnew(iplex,2)=zero
               end do
             end if
!DEBUG_ALTERNATE_ALGO
!          else if (noncoll.and.(lsymnew)) then
!            do mu=1,4
!              do iplex=1,cplex_dij
!                dijnew(iplex,mu)=sumrhoso(iplex,mu)/nsym_used(1)
!                if (abs(dijnew(iplex,mu))<=tol10) dijnew(iplex,mu)=zero
!              end do
!            end do
!DEBUG_ALTERNATE_ALGO
           end if

!          Non-collinear case: store new values of Dij magnetization
           if (noncoll.and.(.not.lsymnew)) then
!            Select on-zero elements
             do mu=1,3
               do iplex=1,cplex_dij
                 rotmag(iplex,mu)=rotmag(iplex,mu)/nsym_used(1)
                 if (abs(rotmag(iplex,mu))<=tol10) rotmag(iplex,mu)=zero
               end do
             end do
!            Transfer back to Dij^{alpha,beta}
             if(.not.lsymnew) then
               dijnew(1,1)=half*(dijnew(1,1)+rotmag(1,3))
               dijnew(2,1)=half*(dijnew(2,1)+rotmag(2,3))
               dijnew(1,2)=      dijnew(1,1)-rotmag(1,3)
               dijnew(2,2)=      dijnew(2,1)-rotmag(2,3)
               dijnew(1,3)=half*(rotmag(1,1)+rotmag(2,2))
               dijnew(2,3)=half*(rotmag(2,1)-rotmag(1,2))
               dijnew(1,4)=half*(rotmag(1,1)-rotmag(2,2))
               dijnew(2,4)=half*(rotmag(2,1)+rotmag(1,2))
             end if
           end if

!          Transfer new value of Dij in suitable pointer
           if (option_dij==0) then
             paw_ij(iatom)%dij(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           else if (option_dij==1) then
             paw_ij(iatom)%dijhat(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           else if (option_dij==2) then
             paw_ij(iatom)%dijU(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           else if (option_dij==3) then
             paw_ij(iatom)%dijxc(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           else if (option_dij==4) then
             paw_ij(iatom)%dijxc_hat(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           else if (option_dij==5) then
             paw_ij(iatom)%dijxc_val(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           else if (option_dij==6) then
             paw_ij(iatom)%dijso(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           else if (option_dij==7) then
             paw_ij(iatom)%dijexxc(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           else if (option_dij==8) then
             paw_ij(iatom)%dijfr(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           else if (option_dij==9) then
             paw_ij(iatom)%dijnd(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           end if

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
     LIBPAW_DEALLOCATE(summag)
     LIBPAW_DEALLOCATE(rotmag)
     LIBPAW_DEALLOCATE(symrec_cart)
     LIBPAW_DEALLOCATE(work)
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

   if (my_natom>0) then
     if(paw_ij(1)%nspden==2.and.paw_ij(1)%nsppol==1) then
       msg='In the antiferromagnetic case, nsym cannot be 1'
       MSG_BUG(msg)
     end if
   end if
   do iatom=1,my_natom
     do ispden=1,paw_ij(iatom)%ndij

       if (option_dij==0) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dij(klmn,ispden))<=tol10) paw_ij(iatom)%dij(klmn,ispden)=zero
         end do
       else if (option_dij==1) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dijhat(klmn,ispden))<=tol10) paw_ij(iatom)%dijhat(klmn,ispden)=zero
         end do
       else if (option_dij==2) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dijU(klmn,ispden))<=tol10) paw_ij(iatom)%dijU(klmn,ispden)=zero
         end do
       else if (option_dij==3) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dijxc(klmn,ispden))<=tol10) paw_ij(iatom)%dijxc(klmn,ispden)=zero
         end do
       else if (option_dij==4) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dijxc_hat(klmn,ispden))<=tol10) paw_ij(iatom)%dijxc_hat(klmn,ispden)=zero
         end do
       else if (option_dij==5) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dijxc_val(klmn,ispden))<=tol10) paw_ij(iatom)%dijxc_val(klmn,ispden)=zero
         end do
       else if (option_dij==6) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dijso(klmn,ispden))<=tol10) paw_ij(iatom)%dijso(klmn,ispden)=zero
         end do
       else if (option_dij==7) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dijexxc(klmn,ispden))<=tol10) paw_ij(iatom)%dijexxc(klmn,ispden)=zero
         end do
       else if (option_dij==8) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dijfr(klmn,ispden))<=tol10) paw_ij(iatom)%dijfr(klmn,ispden)=zero
         end do
       else if (option_dij==9) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dijnd(klmn,ispden))<=tol10) paw_ij(iatom)%dijnd(klmn,ispden)=zero
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
   do iatom=1,my_natom,natinc
     iatom_tot=iatom; if (paral_atom) iatom_tot=my_atmtab(iatom)
     write(msg, '(6a,i3,a)') ch10," PAW TEST:",ch10,&
&     ' ====== Values of ',trim(pertstrg),' in symdij (iatom=',iatom_tot,') (Hartree) ======'
     call wrtout(std_out,msg,wrt_mode)
     optsym=2;if (paw_ij(iatom)%cplex_dij==2.and.ipert>0) optsym=1
     do ispden=1,paw_ij(iatom)%ndij
       if (paw_ij(iatom)%ndij==1) then
         write(msg, '(4a,i3,a)') ch10,&
&         ' *********** ',trim(pertstrg),' (atom ',iatom_tot,') **********'
       else
         write(msg, '(4a,i3,3a)') ch10,&
&         ' *********** ',trim(pertstrg),' (atom ',iatom_tot,', Component ', &
&         trim(dspin(ispden+2*(paw_ij(iatom)%ndij/4))),') **********'
       end if
       call wrtout(std_out,msg,wrt_mode)
       if (paw_ij(iatom)%ndij/=4.or.ispden<=2) then
         call pawio_print_ij(std_out,paw_ij(iatom)%dij(:,ispden),paw_ij(iatom)%lmn2_size,&
&         paw_ij(iatom)%cplex_dij,paw_ij(iatom)%lmn_size,-1,idum,0,pawprtvol,idum,50.d0*dble(3-2*ispden),1,&
&         opt_sym=optsym,mode_paral=wrt_mode)
       else
         if (ipert==0) then
           call pawio_print_ij(std_out,paw_ij(iatom)%dij(:,ispden),paw_ij(iatom)%lmn2_size,&
&           paw_ij(iatom)%cplex_dij,paw_ij(iatom)%lmn_size,-1,idum,0,pawprtvol,idum,50.d0*dble(3-2*ispden),1,&
&           asym_ij=paw_ij(iatom)%dij(:,7-ispden),mode_paral=wrt_mode)
         else
           call pawio_print_ij(std_out,paw_ij(iatom)%dij(:,ispden),paw_ij(iatom)%lmn2_size,&
&           paw_ij(iatom)%cplex_dij,paw_ij(iatom)%lmn_size,-1,idum,0,pawprtvol,idum,50.d0*dble(3-2*ispden),1,&
&           opt_sym=optsym,mode_paral=wrt_mode)
         end if
       end if
     end do
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symdij_symcart'
!End of the abilint section

   implicit none
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
!!  paw_ij(natom)%dij???(cplex_dij*lmn2_size,nspden)=symmetrized dij quantities as output
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symdij_all'
!End of the abilint section

 implicit none

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
 integer,parameter :: MAX_NOPTS=11
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

!FIXME  Dij_hartree and exech_pot are not symmetrized,

 if (ANY(paw_ij(:)%has_dijhartree==2)) then
   msg='symdij does not symmetrize dijhartree term!'
   MSG_WARNING(msg)
   !nopt = nopt + 1
   !options(nopt) = 9
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdij_gather'
!End of the abilint section

 implicit none

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

END MODULE m_pawdij
!!***

!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_denpot
!! NAME
!!  m_paw_denpot
!!
!! FUNCTION
!!  This module contains routines related to PAW on-site densities and on-site potentials.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2018 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_denpot

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_time, only : timab

 use m_pawang,           only : pawang_type
 use m_pawrad,           only : pawrad_type,pawrad_deducer0,poisson,simp_gen
 use m_pawtab,           only : pawtab_type
 use m_paw_an,           only : paw_an_type
 use m_paw_ij,           only : paw_ij_type
 use m_pawfgrtab,        only : pawfgrtab_type
 use m_pawrhoij,         only : pawrhoij_type
 use m_pawdij,           only : pawdijhartree,pawdiju_euijkl,pawdijnd,pawdijso,pawxpot,pawdijfock,symdij,symdij_all
 use m_pawxc,            only : pawxc,pawxc_dfpt,pawxcm,pawxcm_dfpt,pawxcpositron,pawxcmpositron
 use m_paw_finegrid,     only : pawgylm
 use m_paral_atom,       only : get_my_atmtab,free_my_atmtab
 use m_paw_correlations, only : pawuenergy,pawxenergy,setnoccmmp
 use m_paral_atom,       only : get_my_atmtab,free_my_atmtab

 use m_crystal,          only : crystal_t
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype

 implicit none

 private

!public procedures.
 public :: pawdenpot    ! Compute different (PAW) energies, densities and potentials inside PAW spheres
 public :: pawdensities ! Compute PAW on-site densities (all-electron ,pseudo and compensation)
 public :: paw_mknewh0  ! Compute bare PAW on-site Hamiltonian (-> GW calculations)

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_denpot/pawdenpot
!! NAME
!! pawdenpot
!!
!! FUNCTION
!! Compute different (PAW) energies, densities and potentials (or potential-like quantities)
!! inside PAW spheres
!! Can also compute first-order densities potentials and second-order energies (RF calculations).
!!
!! INPUTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  [hyb_mixing, hyb_mixing_sr]= -- optional-- mixing factors for the global (resp. screened) XC hybrid functional
!!  ipert=index of perturbation (used only for RF calculation ; set ipert<=0 for GS calculations.
!!  ixc= choice of exchange-correlation scheme (see above, and below)
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms in unit cell.
!!  nucdipmom(3,my_natom) nuclear dipole moments
!!  nzlmopt= if -1, compute all LM-moments of densities
!!                  initialize "lmselect" (index of non-zero LM-moments of densities)
!!           if  0, compute all LM-moments of densities
!!                  force "lmselect" to .true. (index of non-zero LM-moments of densities)
!!           if  1, compute only non-zero LM-moments of densities (stored before)
!!  option=0: compute both energies and potentials
!!         1: compute only potentials
!!         2: compute only energies
!!  paw_an(my_natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  paw_an0(my_natom) <type(paw_an_type)>=paw arrays given on angular mesh for Ground-State
!                                      used only if ipert>0; must be set equal to paw_an for GS calc.
!!  paw_ij(my_natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawspnorb=flag: 1 if spin-orbit coupling is activated
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  ucvol=unit cell volume (bohr^3)
!!  xclevel= XC functional level
!!  xc_denpos= lowest allowed density (usually for the computation of the XC functionals)
!!  znucl(ntypat)=gives the nuclear charge for all types of atoms
!!
!! OUTPUT
!!  paw_ij(my_natom)%dijhartree(cplex*lmn2_size)=Hartree contribution to dij;
!!                                      Enters into calculation of hartree energy
!!  ==== if option=0 or 2
!!    epaw=contribution to total energy from the PAW "on-site" part
!!    epawdc=contribution to total double-counting energy from the PAW "on-site" part
!!  ==== if option=0 or 2 and ipert<=0
!!    compch_sph=compensation charge integral inside spheres computed over spherical meshes
!!  ==== if (option=0 or 1) and paw_an(:)%has_vxc=1
!!    paw_an(my_natom)%vxc1[m](cplex*mesh_size,:,nspden)=XC potential calculated from "on-site" density
!!    paw_an(my_natom)%vxct1[m](cplex*mesh_size,:,nspden)=XC potential calculated from "on-site" pseudo density
!!    ==== if paw_an(iatom_tot)%has_vxcval==1 compute also XC potentials neglecting core charge
!!      paw_an(my_natom)%vxc1_val[m](cplex*mesh_size,:nspden)=XC potential calculated from spherical valence density
!!      paw_an(my_natom)%vxct1_val[m](cplex*mesh_size,:nspden)=XC potential calculated from spherical valence pseudo density
!!  ==== if nzlmopt==-1,
!!    paw_an(iatom_tot)%lnmselect(lm_size,nspden)=select the non-zero LM-moments of rho1 and trho1
!!  ==== if paw_an(:)%has_vhartree=1
!!    paw_an(my_natom)%vh1(cplex*mesh_size,1,1)=Hartree total potential calculated from "on-site" density
!!  ==== if pawspnorb>0
!!    paw_ij(my_natom)%dijso(cplex*cplex_dij*lmn2_size,nspden)=spin-orbit contribution to dij
!!
!! NOTES
!!  Response function calculations:
!!    In order to compute first- or second-order qunatities, paw_an (resp. paw_ij) datastructures
!!    must contain first-order quantities, namely paw_an1 (resp. paw_ij1).
!!
!! PARENTS
!!      bethe_salpeter,dfpt_scfcv,odamix,paw_qpscgw,respfn,scfcv,screening
!!      sigma
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,pawdensities,pawdijfock,pawdijhartree
!!      pawdijnd,pawdijso,pawrad_deducer0,pawuenergy,pawxc,pawxc_dfpt,pawxcm
!!      pawxcm_dfpt,pawxcmpositron,pawxcpositron,pawxenergy,pawxpot,poisson
!!      setnoccmmp,timab,wrtout,xmpi_sum
!!
!! SOURCE

subroutine pawdenpot(compch_sph,epaw,epawdc,ipert,ixc,&
& my_natom,natom,nspden,ntypat,nucdipmom,nzlmopt,option,paw_an,paw_an0,&
& paw_ij,pawang,pawprtvol,pawrad,pawrhoij,pawspnorb,pawtab,pawxcdev,spnorbscl,xclevel,xc_denpos,ucvol,znucl,&
& electronpositron,mpi_atmtab,comm_atom,vpotzero,hyb_mixing,hyb_mixing_sr) ! optional arguments


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdenpot'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ipert,ixc,my_natom,natom,nspden,ntypat,nzlmopt,option,pawprtvol
 integer,intent(in) :: pawspnorb,pawxcdev,xclevel
 integer,optional,intent(in) :: comm_atom
 real(dp), intent(in) :: spnorbscl,xc_denpos,ucvol
 real(dp),intent(in),optional :: hyb_mixing,hyb_mixing_sr
 real(dp),intent(out) :: compch_sph,epaw,epawdc
 type(electronpositron_type),pointer,optional :: electronpositron
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: nucdipmom(3,my_natom),znucl(ntypat)
 real(dp),intent(out),optional :: vpotzero(2)
 type(paw_an_type),intent(inout) :: paw_an(my_natom)
 type(paw_an_type), intent(in) :: paw_an0(my_natom)
 type(paw_ij_type),intent(inout) :: paw_ij(my_natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: cplex,cplex_dij,cplex_rhoij,has_kxc,has_k3xc,iatom,iatom_tot,idum,ierr,ipositron,irhoij,ispden,itypat,itypat0
 integer :: jrhoij,kklmn,klmn,lm_size,lmn2_size,mesh_size,my_comm_atom,ndij,nkxc1,nk3xc1,nspdiag,nsppol,opt_compch
 integer :: usecore,usepawu,usetcore,usexcnhat,usenhat,usefock
 logical :: cplex_eq_two,keep_vhartree,my_atmtab_allocated,need_kxc,need_k3xc,non_magnetic_xc,paral_atom,pawu_new_algo,temp_vxc
 real(dp) :: e1t10,e1xc,e1xcdc,efock,efockdc,eexc,eexcdc,eexdctemp
 real(dp) :: eexc_val,eexcdc_val,eexex,eexexdc,eextemp,eh2
 real(dp) :: eldaumdc,eldaumdcdc,enucdip,etmp,espnorb,etild1xc,etild1xcdc
 real(dp) :: exccore,exchmix,hyb_mixing_,hyb_mixing_sr_,rdum,ro_im
 character(len=3) :: pertstrg
 character(len=500) :: msg
!arrays
 integer :: idum1(0),idum3(0,0,0)
 integer,pointer :: my_atmtab(:)
 logical,allocatable :: lmselect_cur(:),lmselect_cur_ep(:),lmselect_ep(:),lmselect_tmp(:)
 real(dp) :: ro(2),mpiarr(7),tsec(2)
 real(dp),allocatable :: dij_ep(:),dijfock_vv(:,:),dijfock_cv(:,:),diju_im(:,:)
 real(dp),allocatable :: one_over_rad2(:),kxc_tmp(:,:,:),k3xc_tmp(:,:,:)
 real(dp),allocatable :: nhat1(:,:,:),nhat1_ep(:,:,:)
 real(dp) :: rdum2(0,0),rdum3(0,0,0),rdum3a(0,0,0),rdum4(0,0,0,0)
 real(dp),allocatable :: rho(:),rho1(:,:,:),rho1_ep(:,:,:),rho1xx(:,:,:)
 real(dp),allocatable :: trho1(:,:,:),trho1_ep(:,:,:),vxc_tmp(:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(560,1,tsec)

 if(nzlmopt/=0.and.nzlmopt/=1.and.nzlmopt/=-1) then
   msg='invalid value for variable "nzlmopt".'
   MSG_BUG(msg)
 end if

 usepawu=maxval(pawtab(1:ntypat)%usepawu)
 pawu_new_algo=(usepawu==5.or.usepawu==6)
 if (my_natom>0) then
   if(paw_ij(1)%has_dijhartree==0.and.ipert/=natom+1.and.ipert/=natom+10) then
     msg='dijhartree must be allocated !'
     MSG_BUG(msg)
   end if
   if(paw_ij(1)%has_dijU==0.and.pawu_new_algo.and.ipert/=natom+1.and.ipert/=natom+10) then
     msg='dijU must be allocated !'
     MSG_BUG(msg)
   end if
   if (paw_ij(1)%cplex_rf/=paw_an(1)%cplex) then
     msg='paw_ij()%cplex_rf and paw_an()%cplex must be equal !'
     MSG_BUG(msg)
   end if
   if (pawrhoij(1)%cplex<paw_an(1)%cplex) then
     msg='pawrhoij()%cplex must be >=paw_an()%cplex  !'
     MSG_BUG(msg)
   end if
   if (ipert<=0.and.paw_ij(1)%cplex_rf/=1) then
     msg='paw_ij()%cplex_rf must be 1 for GS calculations !'
     MSG_BUG(msg)
   end if
   if (ipert>0.and.(ipert<=natom.or.ipert==natom+2).and.paw_an0(1)%has_kxc/=2) then
     msg='XC kernels for ground state must be in memory !'
     MSG_BUG(msg)
   end if
   if(paw_an(1)%has_vxc==0.and.(option==0.or.option==1).and. &
&   .not.(ipert==natom+1.or.ipert==natom+10)) then
     msg='vxc1 and vxct1 must be allocated !'
     MSG_BUG(msg)
   end if
   if (ipert>0.and.paw_an(1)%has_vhartree==1) then
     msg='computation of vhartree not compatible with RF (ipert>0) !'
     MSG_BUG(msg)
   end if
   if (ipert>0.and.paw_an(1)%has_vxcval==1.and.(option==0.or.option==1)) then
     msg='computation of vxc_val not compatible with RF (ipert>0) !'
     MSG_BUG(msg)
   end if
 end if

 ipositron=0
 if (present(electronpositron)) then
   ipositron=electronpositron_calctype(electronpositron)
   if (ipositron==1.and.pawtab(1)%has_kij/=2) then
     msg='kij must be in memory for electronpositron%calctype=1 !'
     MSG_BUG(msg)
   end if
   if (ipert>0) then
     msg='electron-positron calculation not available for ipert>0 !'
     MSG_ERROR(msg)
   end if
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,&
& my_natom_ref=my_natom)

!For some perturbations, nothing to do
 if (ipert==natom+1.or.ipert==natom+10) then
   if (option/=1) then
     epaw=zero;epawdc=zero
   end if
   return
 end if

!Various inits
 hyb_mixing_   =zero ; if(present(hyb_mixing))    hyb_mixing_   =hyb_mixing
 hyb_mixing_sr_=zero ; if(present(hyb_mixing_sr)) hyb_mixing_sr_=hyb_mixing_sr
 usefock=0;if (abs(hyb_mixing_)>tol8.or.abs(hyb_mixing_sr_)>tol8) usefock=1
 usexcnhat=maxval(pawtab(1:ntypat)%usexcnhat)
 usenhat = usexcnhat
 keep_vhartree=(maxval(paw_an(:)%has_vhartree)>0)
 if(keep_vhartree) usenhat = 1
 non_magnetic_xc=(usepawu==4).or.(usepawu==14)
 compch_sph=-1.d5
 opt_compch=0;if (option/=1.and.ipert<=0) opt_compch=1
 if (opt_compch==1) compch_sph=zero
 nspdiag=1;if (nspden==2) nspdiag=2
 nsppol=1;if (my_natom>0) nsppol=pawrhoij(1)%nsppol
 pertstrg=" ";if (ipert>0) pertstrg="(1)"
 ro(:)=zero

!Init energies
 if (option/=1) then
   e1xc=zero     ; e1xcdc=zero
   etild1xc=zero ; etild1xcdc=zero
   exccore=zero  ; eh2=zero ; e1t10=zero
   eldaumdc=zero ; eldaumdcdc=zero
   eexex=zero    ; eexexdc=zero
   eextemp=zero  ; eexdctemp=zero
   espnorb=zero  ; enucdip=zero
   efock=zero    ; efockdc=zero
   if (ipositron/=0) then
     electronpositron%e_paw  =zero
     electronpositron%e_pawdc=zero
   end if
 end if
!vpotzero is needed for both the energy and the potential
 if (present(vpotzero)) vpotzero(:)=zero

!if PAW+U, compute noccmmp^{\sigma}_{m,m'} occupation matrix
 if (usepawu>0.and.ipert<=0.and.ipositron/=1) then
   if (paral_atom) then
     call setnoccmmp(1,0,rdum4,0,0,idum3,my_natom,natom,0,1,nsppol,0,ntypat,&
&     paw_ij,pawang,pawprtvol,pawrhoij,pawtab,rdum2,idum1,idum1,0,usepawu,&
&     comm_atom=my_comm_atom,mpi_atmtab=mpi_atmtab)
   else
     call setnoccmmp(1,0,rdum4,0,0,idum3,my_natom,natom,0,1,nsppol,0,ntypat,&
&     paw_ij,pawang,pawprtvol,pawrhoij,pawtab,rdum2,idum1,idum1,0,usepawu)
   end if
 end if

!Print some titles
 if (abs(pawprtvol)>=2) then
   if (nzlmopt<1) write(msg, '(6a)') ch10,' PAW TEST:',ch10,&
   ' ====== Moments of (n1-tn1)',trim(pertstrg),' ========='
   if (nzlmopt==1) write(msg, '(6a)') ch10,' PAW TEST:',ch10,&
   ' ==== Non-zero Moments of (n1-tn1)',trim(pertstrg),' ===='
   call wrtout(std_out,msg,'COLL')
   if (usexcnhat/=0) then
     write(msg, '(6a)')' The moments of (n1-tn1-nhat1)',trim(pertstrg),' must be very small...'
     call wrtout(std_out,msg,'COLL')
   end if
 end if

!================ Big loop on atoms =======================
!==========================================================

 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
   itypat=pawrhoij(iatom)%itypat
   exchmix=pawtab(itypat)%exchmix
   lmn2_size=paw_ij(iatom)%lmn2_size
   lm_size=paw_an(iatom)%lm_size
   mesh_size=pawtab(itypat)%mesh_size

   usecore=1;usetcore =pawtab(itypat)%usetcore
   if (ipert/=0) usecore=0  ! This is true for phonons and Efield pert.
   if (ipert/=0) usetcore=0 ! This is true for phonons and Efield pert.
   has_kxc =paw_an(iatom)%has_kxc ;need_kxc =(has_kxc ==1)
   has_k3xc=paw_an(iatom)%has_k3xc;need_k3xc=(has_k3xc==1)
   cplex=paw_an(iatom)%cplex
   cplex_dij=paw_ij(iatom)%cplex_dij
   cplex_rhoij=pawrhoij(iatom)%cplex
   ndij=paw_ij(iatom)%ndij

!  Allocations of "on-site" densities
   ABI_ALLOCATE(rho1 ,(cplex*mesh_size,lm_size,nspden))
   ABI_ALLOCATE(trho1,(cplex*mesh_size,lm_size,nspden))
   ABI_ALLOCATE(nhat1,(cplex*mesh_size,lm_size,nspden*usenhat))
   rho1(:,:,:)=zero;trho1(:,:,:)=zero;nhat1(:,:,:)=zero
   if (ipositron/=0) then ! Additional allocation for the electron-positron case
     ABI_ALLOCATE(rho1_ep ,(cplex*mesh_size,lm_size,nspden))
     ABI_ALLOCATE(trho1_ep,(cplex*mesh_size,lm_size,nspden))
     ABI_ALLOCATE(nhat1_ep,(cplex*mesh_size,lm_size,nspden*usenhat))
   end if
   ABI_ALLOCATE(lmselect_cur,(lm_size))
   lmselect_cur(:)=.true.
   if (nzlmopt==1) lmselect_cur(:)=paw_an(iatom)%lmselect(:)

!  Store some usefull quantities
   itypat0=0;if (iatom>1) itypat0=pawrhoij(iatom-1)%itypat
   if (itypat/=itypat0) then
     ABI_ALLOCATE(one_over_rad2,(mesh_size))
     one_over_rad2(1)=zero
     one_over_rad2(2:mesh_size)=one/pawrad(itypat)%rad(2:mesh_size)**2
   end if

!  Need to allocate vxc1 in particular cases
   if (pawspnorb>0.and.ipert==0.and.option==2.and.ipositron/=1.and. &
&      cplex_rhoij==2.and.paw_an(iatom)%has_vxc==0) then
!    These should already be allocated in paw_an_init!
     if (allocated(paw_an(iatom)%vxc1))  then
       ABI_DEALLOCATE(paw_an(iatom)%vxc1)
     end if
     if (pawxcdev==0)then
       ABI_ALLOCATE(paw_an(iatom)%vxc1,(cplex*mesh_size,paw_an(iatom)%angl_size,nspden))
     else
       ABI_ALLOCATE(paw_an(iatom)%vxc1,(cplex*mesh_size,lm_size,nspden))
     end if
     paw_an(iatom)%has_vxc=1
     temp_vxc=.true.
   else
     temp_vxc=.false.
   end if

!  ===== Compute "on-site" densities (n1, ntild1, nhat1) =====
!  ==========================================================

   call pawdensities(compch_sph,cplex,iatom_tot,lmselect_cur,paw_an(iatom)%lmselect,lm_size,&
&   nhat1,nspden,nzlmopt,opt_compch,1-usenhat,-1,1,pawang,pawprtvol,pawrad(itypat),&
&   pawrhoij(iatom),pawtab(itypat),rho1,trho1,one_over_rad2=one_over_rad2)

   if (ipositron/=0) then
!    Electron-positron calculation: need additional on-site densities:
!    if ipositron==1, need electronic on-site densities
!    if ipositron==2, need positronic on-site densities
     ABI_ALLOCATE(lmselect_ep,(lm_size))
     ABI_ALLOCATE(lmselect_cur_ep,(lm_size))
     lmselect_cur_ep(:)=.true.
     if (nzlmopt==1) lmselect_cur_ep(:)=electronpositron%lmselect_ep(1:lm_size,iatom)

     call pawdensities(rdum,cplex,iatom_tot,lmselect_cur_ep,lmselect_ep,&
&     lm_size,nhat1_ep,nspden,nzlmopt,0,1-usenhat,-1,0,pawang,0,pawrad(itypat),&
&     electronpositron%pawrhoij_ep(iatom),pawtab(itypat),&
&     rho1_ep,trho1_ep,one_over_rad2=one_over_rad2)

     if (nzlmopt<1) electronpositron%lmselect_ep(1:lm_size,iatom)=lmselect_ep(1:lm_size)
     ABI_DEALLOCATE(lmselect_ep)
     ABI_DEALLOCATE(lmselect_cur_ep)
   end if

!  =========== Compute XC potentials and energies ===========
!  ==========================================================

!  Temporary storage
   nkxc1 =0;if (paw_an(iatom)%has_kxc /=0) nkxc1 =paw_an(iatom)%nkxc1
   nk3xc1=0;if (paw_an(iatom)%has_k3xc/=0.and.pawxcdev==0) nk3xc1=paw_an(iatom)%nk3xc1
   if (pawxcdev/=0) then
     ABI_ALLOCATE(vxc_tmp,(cplex*mesh_size,lm_size,nspden))
     if (need_kxc) then
       ABI_ALLOCATE(kxc_tmp,(mesh_size,lm_size,nkxc1))
     end if
     if (need_k3xc) then
       msg = 'Computation of k3xc with pawxcdev/=0 is not implemented yet!'
       MSG_BUG(msg)
     end if
   end if
   if (pawxcdev==0) then
     ABI_ALLOCATE(vxc_tmp,(cplex*mesh_size,pawang%angl_size,nspden))
     vxc_tmp(:,:,:)=zero
     if (need_kxc) then
       ABI_ALLOCATE(kxc_tmp,(mesh_size,pawang%angl_size,nkxc1))
     end if
     if (need_k3xc) then
       ABI_ALLOCATE(k3xc_tmp,(mesh_size,pawang%angl_size,nk3xc1))
     end if
   end if
   idum=0
   if (.not.allocated(kxc_tmp))  then
     ABI_ALLOCATE(kxc_tmp,(0,0,0))
   end if
   if (.not.allocated(k3xc_tmp))  then
     ABI_ALLOCATE(k3xc_tmp,(0,0,0))
   end if

!  ===== Vxc1 term =====
   if (ipositron/=1) then
     if (pawxcdev/=0) then
       if (ipert==0) then
         call pawxcm(pawtab(itypat)%coredens,eexc,eexcdc,idum,ixc,kxc_tmp,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,nkxc1,non_magnetic_xc,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),pawxcdev,rho1,usecore,0,vxc_tmp,xclevel,xc_denpos)
       else
         call pawxcm_dfpt(pawtab(itypat)%coredens,cplex,cplex,eexc,ixc,paw_an0(iatom)%kxc1,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,paw_an0(iatom)%nkxc1,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),rho1,usecore,0,vxc_tmp,xclevel)
         eexcdc=zero
       end if
     else
       if (ipert==0) then
         call pawxc(pawtab(itypat)%coredens,eexc,eexcdc,ixc,kxc_tmp,k3xc_tmp,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,nkxc1,nk3xc1,non_magnetic_xc,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),rho1,usecore,0,vxc_tmp,xclevel,xc_denpos)
       else
         call pawxc_dfpt(pawtab(itypat)%coredens,cplex,cplex,eexc,ixc,paw_an0(iatom)%kxc1,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,paw_an0(iatom)%nkxc1,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),rho1,usecore,0,paw_an0(iatom)%vxc1,vxc_tmp,xclevel)
         eexcdc=zero
       end if
     end if

     if (option/=1) then
       e1xc=e1xc+eexc
       e1xcdc=e1xcdc+eexcdc
     end if
     if (option<2.or.temp_vxc) paw_an(iatom)%vxc1(:,:,:)=vxc_tmp(:,:,:)
     if (need_kxc .and.nkxc1>0 ) paw_an(iatom)%kxc1(:,:,:) =kxc_tmp(:,:,:)
     if (need_k3xc.and.nk3xc1>0) paw_an(iatom)%k3xc1(:,:,:)=k3xc_tmp(:,:,:)
   else ! ipositron==1
     if (option<2.or.temp_vxc) paw_an(iatom)%vxc1(:,:,:)=zero
     if (need_kxc.and.nkxc1>0) paw_an(iatom)%kxc1(:,:,:)=zero
   end if


!  Additional electron-positron XC term (if ipositron/=0)
   if (ipositron/=0) then
     if (pawxcdev/=0) then
       call pawxcmpositron(ipositron,pawtab(itypat)%coredens,eexc,eexcdc,electronpositron%ixcpositron,&
&       lm_size,paw_an(iatom)%lmselect,electronpositron%lmselect_ep(1:lm_size,iatom),&
&       nhat1,nhat1_ep,mesh_size,nspden,option,pawang,pawrad(itypat),pawxcdev,&
&       electronpositron%posdensity0_limit,rho1,rho1_ep,usecore,0,vxc_tmp,xc_denpos)
     else
       call pawxcpositron(ipositron,pawtab(itypat)%coredens,eexc,eexcdc,electronpositron%ixcpositron,&
&       lm_size,paw_an(iatom)%lmselect,electronpositron%lmselect_ep(1:lm_size,iatom),&
&       nhat1,nhat1_ep,mesh_size,nspden,option,pawang,pawrad(itypat),&
&       electronpositron%posdensity0_limit,rho1,rho1_ep,usecore,0,vxc_tmp,xc_denpos)
     end if
     if (option/=1) then
       electronpositron%e_paw  =electronpositron%e_paw  +eexc
       electronpositron%e_pawdc=electronpositron%e_pawdc+eexcdc
     end if
     if (option<2.or.temp_vxc) paw_an(iatom)%vxc1(:,:,:)=paw_an(iatom)%vxc1(:,:,:)+vxc_tmp(:,:,:)
     if (need_kxc.and.nkxc1>0) paw_an(iatom)%kxc1(:,:,:)=paw_an(iatom)%kxc1(:,:,:)+kxc_tmp(:,:,:)
   end if

!  ===== tVxc1 term =====
   if (ipositron/=1) then
     if (pawxcdev/=0) then
       if (ipert==0) then
         call pawxcm(pawtab(itypat)%tcoredens(:,1),&
&         eexc,eexcdc,idum,ixc,kxc_tmp,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,nkxc1,non_magnetic_xc,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),pawxcdev,trho1,usetcore,2*usexcnhat,vxc_tmp,xclevel,xc_denpos)
       else
         call pawxcm_dfpt(pawtab(itypat)%tcoredens(:,1),&
&         cplex,cplex,eexc,ixc,paw_an0(iatom)%kxct1,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,paw_an0(iatom)%nkxc1,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),trho1,usetcore,2*usexcnhat,vxc_tmp,xclevel)
         eexcdc=zero
       end if
     else
       if (ipert==0) then
         call pawxc(pawtab(itypat)%tcoredens(:,1),&
&         eexc,eexcdc,ixc,kxc_tmp,k3xc_tmp,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,nkxc1,nk3xc1,non_magnetic_xc,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),trho1,usetcore,2*usexcnhat,vxc_tmp,xclevel,xc_denpos)
       else
         call pawxc_dfpt(pawtab(itypat)%tcoredens(:,1),&
&         cplex,cplex,eexc,ixc,paw_an0(iatom)%kxct1,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,paw_an0(iatom)%nkxc1,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),trho1,usetcore,2*usexcnhat,paw_an0(iatom)%vxct1,vxc_tmp,xclevel)
         eexcdc=zero
       end if
     end if
     if (option/=1) then
       etild1xc=etild1xc+eexc
       etild1xcdc=etild1xcdc+eexcdc
     end if
     if (option<2) paw_an(iatom)%vxct1(:,:,:)=vxc_tmp(:,:,:)
     if (need_kxc.and. nkxc1>0 ) paw_an(iatom)%kxct1(:,:,:) =kxc_tmp(:,:,:)
     if (need_k3xc.and.nk3xc1>0) paw_an(iatom)%k3xct1(:,:,:)=k3xc_tmp(:,:,:)
   else ! ipositron==1
     if (option<2) paw_an(iatom)%vxct1(:,:,:)=zero
     if (need_kxc.and.nkxc1>0) paw_an(iatom)%kxct1(:,:,:)=zero
   end if

!  Additional electron-positron XC term (if ipositron/=0)
   if (ipositron/=0) then
     if (pawxcdev/=0) then
       call pawxcmpositron(ipositron,pawtab(itypat)%tcoredens(:,1),&
&       eexc,eexcdc,electronpositron%ixcpositron,&
&       lm_size,paw_an(iatom)%lmselect,electronpositron%lmselect_ep(1:lm_size,iatom),&
&       nhat1,nhat1_ep,mesh_size,nspden,option,pawang,pawrad(itypat),pawxcdev,&
&       electronpositron%posdensity0_limit,trho1,trho1_ep,usetcore,2*usexcnhat,vxc_tmp,xc_denpos)
     else
       call pawxcpositron(ipositron,pawtab(itypat)%tcoredens(:,1),&
&       eexc,eexcdc,electronpositron%ixcpositron,&
&       lm_size,paw_an(iatom)%lmselect,electronpositron%lmselect_ep(1:lm_size,iatom),&
&       nhat1,nhat1_ep,mesh_size,nspden,option,pawang,pawrad(itypat),&
&       electronpositron%posdensity0_limit,trho1,trho1_ep,usetcore,2*usexcnhat,vxc_tmp,xc_denpos)
     end if
     if (option/=1) then
       electronpositron%e_paw  =electronpositron%e_paw  -eexc
       electronpositron%e_pawdc=electronpositron%e_pawdc-eexcdc
     end if
     if (option<2) paw_an(iatom)%vxct1(:,:,:)=paw_an(iatom)%vxct1(:,:,:)+vxc_tmp(:,:,:)
     if (need_kxc.and.nkxc1>0) paw_an(iatom)%kxct1(:,:,:)=paw_an(iatom)%kxct1(:,:,:)+kxc_tmp(:,:,:)
   end if

!  Update flags defining the state of vxc and kxc
   if (option<2) paw_an(iatom)%has_vxc=2
   if (need_kxc.and.nkxc1>0) paw_an(iatom)%has_kxc=2

!  Update core XC conjtribution to energy
   if (option/=1.and.ipositron/=1) exccore=exccore+pawtab(itypat)%exccore
   if (ipositron/=0)  then
     ABI_DEALLOCATE(rho1_ep)
     ABI_DEALLOCATE(trho1_ep)
     ABI_DEALLOCATE(nhat1_ep)
   end if

!  =========== Compute valence-only XC potentials ===========
!  ==========================================================
   if (ipert==0.and.paw_an(iatom)%has_vxcval==1.and.(option==0.or.option==1)) then
     if (.not.allocated(paw_an(iatom)%vxc1_val).or..not.allocated(paw_an(iatom)%vxct1_val)) then
       msg=' vxc1_val and vxct1_val must be associated'
       MSG_BUG(msg)
     end if
!    ===== Vxc1_val term, vxc[n1] =====
     if (pawxcdev/=0) then
       write(msg,'(4a,es16.6)')ch10,&
&       ' pawdenpot : Computing valence-only v_xc[n1] using moments ',ch10,&
&       '             Min density rho1 = ',MINVAL(rho1)
       call wrtout(std_out,msg,'COLL')
       call pawxcm(pawtab(itypat)%coredens,eexc_val,eexcdc_val,idum,ixc,kxc_tmp,lm_size,&
&       paw_an(iatom)%lmselect,nhat1,nkxc1,non_magnetic_xc,mesh_size,nspden,option,&
&       pawang,pawrad(itypat),pawxcdev,rho1,0,0,vxc_tmp,xclevel,xc_denpos)
     else
       write(msg,'(2a)')ch10,' pawdenpot : Computing valence-only v_xc[n1] using angular mesh '
       call wrtout(std_out,msg,'COLL')

       call pawxc(pawtab(itypat)%coredens,eexc_val,eexcdc_val,ixc,kxc_tmp,k3xc_tmp,lm_size,&
&       paw_an(iatom)%lmselect,nhat1,nkxc1,nk3xc1,non_magnetic_xc,mesh_size,nspden,option,&
&       pawang,pawrad(itypat),rho1,0,0,vxc_tmp,xclevel,xc_denpos)
     end if
     if (option<2) paw_an(iatom)%vxc1_val(:,:,:)=vxc_tmp(:,:,:)

!    ===== tVxc1_val term =====
     if (pawxcdev/=0) then
       if (usexcnhat/=0) then
         write(msg,'(4a,e16.6,2a,es16.6)')ch10,&
&         ' pawdenpot : Computing valence-only v_xc[tn1+nhat] using moments ',ch10,&
&         '             Min density trho1        = ',MINVAL(trho1),ch10,&
&         '             Min density trho1 + nhat = ',MINVAL(trho1+nhat1)
       else
         write(msg,'(4a,e16.6)')ch10,&
&         ' pawdenpot : Computing valence-only v_xc[tn1] using moments ',ch10,&
&         '             Min density trho1        = ',MINVAL(trho1)
       end if
       call wrtout(std_out,msg,'COLL')
       call pawxcm(pawtab(itypat)%tcoredens(:,1),&
&       eexc_val,eexcdc_val,idum,ixc,kxc_tmp,lm_size,&
&       paw_an(iatom)%lmselect,nhat1,nkxc1,non_magnetic_xc,mesh_size,nspden,option,&
&       pawang,pawrad(itypat),pawxcdev,trho1,0,2*usexcnhat,vxc_tmp,xclevel,xc_denpos)
     else
       write(msg,'(2a)')ch10,' pawdenpot : Computing valence-only v_xc[tn1+nhat] using angular mesh'
       call wrtout(std_out,msg,'COLL')
       call pawxc(pawtab(itypat)%tcoredens(:,1),&
&       eexc_val,eexcdc_val,ixc,kxc_tmp,k3xc_tmp,lm_size,&
&       paw_an(iatom)%lmselect,nhat1,nkxc1,nk3xc1,non_magnetic_xc,mesh_size,nspden,option,&
&       pawang,pawrad(itypat),trho1,0,2*usexcnhat,vxc_tmp,xclevel,xc_denpos)
     end if
     if (option<2) then
       paw_an(iatom)%vxct1_val(:,:,:)=vxc_tmp(:,:,:)
       paw_an(iatom)%has_vxcval=2
     end if
   end if ! valence-only XC potentials

   ABI_DEALLOCATE(vxc_tmp)
   ABI_DEALLOCATE(kxc_tmp)
   ABI_DEALLOCATE(k3xc_tmp)

!  ===== Compute first part of local exact-exchange energy term =====
!  ===== Also compute corresponding potential                   =====
!  ==================================================================

   if (pawtab(itypat)%useexexch>0.and.ipert==0.and.ipositron/=1) then

!    ===== Re-compute a partial "on-site" density n1 (only l=lexexch contrib.)
     ABI_ALLOCATE(rho1xx,(mesh_size,lm_size,nspden))
     ABI_ALLOCATE(lmselect_tmp,(lm_size))
     lmselect_tmp(:)=lmselect_cur(:)
     call pawdensities(rdum,cplex,iatom_tot,lmselect_cur,lmselect_tmp,lm_size,rdum3,nspden,&
&     1,0,2,pawtab(itypat)%lexexch,0,pawang,pawprtvol,pawrad(itypat),&
&     pawrhoij(iatom),pawtab(itypat),rho1xx,rdum3a,one_over_rad2=one_over_rad2)
     ABI_DEALLOCATE(lmselect_tmp)
!    ===== Re-compute Exc1 and Vxc1; for local exact-exchange, this is done in GGA only
     ABI_ALLOCATE(vxc_tmp,(mesh_size,lm_size,nspden))
     ABI_ALLOCATE(kxc_tmp,(mesh_size,lm_size,nkxc1))
     call pawxcm(pawtab(itypat)%coredens,eextemp,eexdctemp,pawtab(itypat)%useexexch,ixc,kxc_tmp,lm_size,&
&     paw_an(iatom)%lmselect,nhat1,nkxc1,non_magnetic_xc,mesh_size,nspden,option,pawang,pawrad(itypat),pawxcdev,&
&     rho1xx,0,0,vxc_tmp,xclevel,xc_denpos)
     if (option/=1) then
       e1xc=e1xc-eextemp*exchmix
       e1xcdc=e1xcdc-eexdctemp*exchmix
     end if
     if (option<2) then
       paw_an(iatom)%vxc_ex(:,:,:)=vxc_tmp(:,:,:)
       paw_an(iatom)%has_vxc_ex=2
     end if
     ABI_DEALLOCATE(rho1xx)
     ABI_DEALLOCATE(vxc_tmp)
     ABI_DEALLOCATE(kxc_tmp)

   end if ! useexexch

   itypat0=0;if (iatom<my_natom) itypat0=pawrhoij(iatom+1)%itypat
   if (itypat/=itypat0) then
     ABI_DEALLOCATE(one_over_rad2)
   end if

   ABI_DEALLOCATE(lmselect_cur)

!  ==== Compute Hartree potential terms and some energy terms ====
!  ===============================================================

!  Electron-positron calculation: compute Dij due to fixed particles (elec. or pos. depending on calctype)
   if (ipositron/=0) then
     ABI_CHECK(cplex==1,'BUG in pawdenpot: cplex should be 1 for electron-positron!')
     ABI_ALLOCATE(dij_ep,(cplex*lmn2_size))
     call pawdijhartree(cplex,dij_ep,nspden,electronpositron%pawrhoij_ep(iatom),pawtab(itypat))
     if (option/=1) then
       do ispden=1,nspdiag
         jrhoij=1
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=pawrhoij(iatom)%rhoijselect(irhoij)
           kklmn=klmn;if (cplex==2) kklmn=2*klmn-1
           ro(1)=pawrhoij(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
           electronpositron%e_paw  =electronpositron%e_paw  -ro(1)*dij_ep(kklmn)
           electronpositron%e_pawdc=electronpositron%e_pawdc-ro(1)*dij_ep(kklmn)
           if (ipositron==1) e1t10=e1t10+ro(1)*two*(pawtab(itypat)%kij(klmn)-pawtab(itypat)%dij0(klmn))
           jrhoij=jrhoij+pawrhoij(iatom)%cplex
         end do
       end do
     end if
   end if

!  Hartree Dij computation
   if (ipositron/=1) then
     call pawdijhartree(cplex,paw_ij(iatom)%dijhartree,nspden,pawrhoij(iatom),pawtab(itypat))
   else
     paw_ij(iatom)%dijhartree(:)=zero
   end if
   paw_ij(iatom)%has_dijhartree=2

!  Hartree energy computation
   if (option/=1) then
     if (cplex==1.or.ipert==0) then
       do ispden=1,nspdiag
         jrhoij=1
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=pawrhoij(iatom)%rhoijselect(irhoij)
           ro(1)=pawrhoij(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
           eh2=eh2    +ro(1)*paw_ij(iatom)%dijhartree(klmn)
           e1t10=e1t10+ro(1)*pawtab(itypat)%dij0(klmn)
           jrhoij=jrhoij+pawrhoij(iatom)%cplex
         end do
       end do
     else
       do ispden=1,nspdiag
         jrhoij=1
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=pawrhoij(iatom)%rhoijselect(irhoij);kklmn=klmn+lmn2_size
           ro(1:2)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+1,ispden)*pawtab(itypat)%dltij(klmn)
           eh2=eh2+ro(1)*paw_ij(iatom)%dijhartree(klmn)+ro(2)*paw_ij(iatom)%dijhartree(kklmn)
!          Imaginary part (not used)
!          eh2=eh2+ro(2)*paw_ij(iatom)%dijhartree(klmn)-ro(1)*paw_ij(iatom)%dijhartree(kklmn)
           jrhoij=jrhoij+pawrhoij(iatom)%cplex
         end do
       end do
     end if
   end if

!  Electron-positron calculation: add electron and positron
   if (ipositron/=0) then
     paw_ij(iatom)%dijhartree(:)=paw_ij(iatom)%dijhartree(:)-dij_ep(:)
     ABI_DEALLOCATE(dij_ep)
   end if

!  Compute 1st moment of total Hartree potential VH(n_Z+n_core+n1)
!  equation 10 (density) and up to 43 (Hartree potential of density)
!    of Kresse and Joubert PRB 59 1758 (1999) [[cite:Kresse1999]]
   keep_vhartree=(paw_an(iatom)%has_vhartree>0)
   if ((pawspnorb>0.and.ipert==0.and.ipositron/=1).or.keep_vhartree) then

     ABI_CHECK(cplex==1,'BUG in pawdenpot: vh1 not available for cplex=2!')

     !In the first clause case, would it not be simpler just to turn on has_vhartree?
     if (.not. allocated(paw_an(iatom)%vh1)) then
       ABI_ALLOCATE(paw_an(iatom)%vh1,(cplex*mesh_size,1,1))
     end if
     if (.not. allocated(paw_an(iatom)%vht1)) then
       ABI_ALLOCATE(paw_an(iatom)%vht1,(cplex*mesh_size,1,1))
     end if

     !Construct vh1
     !  The sqrt(4pi) factor comes from the fact we are calculating the spherical moments,
     !   and for the 00 channel the prefactor of Y_00 = 2 sqrt(pi)
     ABI_ALLOCATE(rho,(mesh_size))
     rho(1:mesh_size)=(rho1(1:mesh_size,1,1) + sqrt(four_pi)*pawtab(itypat)%coredens(1:mesh_size)) &
&                    *four_pi*pawrad(itypat)%rad(1:mesh_size)**2
!    rho(1:mesh_size)=rho1(1:mesh_size,1,1)*four_pi*pawrad(itypat)%rad(1:mesh_size)**2
     call poisson(rho,0,pawrad(itypat),paw_an(iatom)%vh1(:,1,1))
     paw_an(iatom)%vh1(2:mesh_size,1,1)=(paw_an(iatom)%vh1(2:mesh_size,1,1) &
&                                      -sqrt(four_pi)*znucl(itypat))/pawrad(itypat)%rad(2:mesh_size)
!    TODO: check this is equivalent to the previous version (commented) which explicitly recalculated VH(coredens)
!    DONE: numerically there are residual differences on abiref (7th digit).
!         paw_an(iatom)%vh1(2:mesh_size,1,1)=paw_an(iatom)%vh1(2:mesh_size,1,1)/pawrad(itypat)%rad(2:mesh_size) &
!&                                          +sqrt(four_pi) * pawtab(itypat)%VHnZC(2:mesh_size)
     call pawrad_deducer0(paw_an(iatom)%vh1(:,1,1),mesh_size,pawrad(itypat))

!    !Same for vht1
     rho = zero
     if (usenhat /= 0) rho(1:mesh_size)=nhat1(1:mesh_size,1,1)
     rho(1:mesh_size)=(rho(1:mesh_size) + trho1(1:mesh_size,1,1) + sqrt(four_pi)*pawtab(itypat)%tcoredens(1:mesh_size,1)) &
&                    *four_pi*pawrad(itypat)%rad(1:mesh_size)**2
!    rho(1:mesh_size)=(rho(1:mesh_size) + trho1(1:mesh_size,1,1))*four_pi*pawrad(itypat)%rad(1:mesh_size)**2
     call poisson(rho,0,pawrad(itypat),paw_an(iatom)%vht1(:,1,1))
     paw_an(iatom)%vht1(2:mesh_size,1,1)=(paw_an(iatom)%vht1(2:mesh_size,1,1) &
&     -sqrt(four_pi)*znucl(itypat))/pawrad(itypat)%rad(2:mesh_size)
!    paw_an(iatom)%vht1(2:mesh_size,1,1)=paw_an(iatom)%vht1(2:mesh_size,1,1)/pawrad(itypat)%rad(2:mesh_size) &
!&                                      +sqrt(four_pi) * pawtab(itypat)%vhtnzc(2:mesh_size)
     call pawrad_deducer0(paw_an(iatom)%vht1(:,1,1),mesh_size,pawrad(itypat))

     paw_an(iatom)%has_vhartree=2
     ABI_DEALLOCATE(rho)
   end if

!  ========= Compute PAW+U and energy contribution  =========
!  ==========================================================

   if (pawtab(itypat)%usepawu>0.and.ipositron/=1.and.option/=1.and.pawtab(itypat)%usepawu<10) then
     if (ipert==0.and.option/=1.and.(.not.pawu_new_algo)) then
       call pawuenergy(iatom_tot,eldaumdc,eldaumdcdc,paw_ij(iatom)%noccmmp, &
&                      paw_ij(iatom)%nocctot,pawprtvol,pawtab(itypat))
     else
!      PAW+U Dij computation from euijkl
       ABI_CHECK(cplex_dij==1,'cplex_dij not yet available!')
       ABI_CHECK(ndij/=4,'ndij not yet available!')
       cplex_eq_two=.true.; if (paw_ij(iatom)%cplex_rf==1.or.ipert==0) cplex_eq_two=.false.
       if (cplex_eq_two.and.option/=1) then
         ABI_ALLOCATE(diju_im,(pawtab(itypat)%lmn2_size,ndij))
         call pawdiju_euijkl(cplex,cplex_dij,paw_ij(iatom)%dijU,ndij, &
&                            pawrhoij(iatom),pawtab(itypat),diju_im=diju_im)
       else
         call pawdiju_euijkl(cplex,cplex_dij,paw_ij(iatom)%dijU,ndij, &
&                            pawrhoij(iatom),pawtab(itypat))
       end if
       paw_ij(iatom)%has_dijU=2
!      PAW+U energy computation
       if (option/=1) then
         if (.not.cplex_eq_two) then
           do ispden=1,ndij
             jrhoij=1
             do irhoij=1,pawrhoij(iatom)%nrhoijsel
               klmn=pawrhoij(iatom)%rhoijselect(irhoij)
               ro(1)=pawrhoij(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
               etmp = half*ro(1)*paw_ij(iatom)%dijU(klmn,ispden)
               eldaumdc   = eldaumdc   + etmp
               eldaumdcdc = eldaumdcdc - etmp
               jrhoij=jrhoij+cplex_rhoij
             end do
           end do
         else
           do ispden=1,ndij
             jrhoij=1
             do irhoij=1,pawrhoij(iatom)%nrhoijsel
               klmn=pawrhoij(iatom)%rhoijselect(irhoij);kklmn=2*klmn-1
               ro(1:2)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+1,ispden)*pawtab(itypat)%dltij(klmn)
               ro_im = pawrhoij(iatom)%rhoijim(klmn,ispden)
               if (abs(pawtab(itypat)%dltij(klmn)-1)<tol14) then ! case i==j
                 etmp = half*(ro(1)*paw_ij(iatom)%dijU(kklmn,ispden)+(ro(2)+ro_im)*paw_ij(iatom)%dijU(kklmn+1,ispden))
               else ! case i/=j (we remove the factor half to take into account that dltij=2)
                 etmp =        ro(1)*paw_ij(iatom)%dijU(kklmn  ,ispden)                        ! re  * re
                 etmp = etmp + ro_im*diju_im(klmn,ispden)                                      ! ima * ima
                 etmp = etmp + ro(2)*(paw_ij(iatom)%dijU(kklmn+1,ispden)-diju_im(klmn,ispden)) ! imb * imb
               end if
               eldaumdc   = eldaumdc   + etmp
               eldaumdcdc = eldaumdcdc - etmp
               jrhoij=jrhoij+cplex_rhoij
             end do
           end do
           ABI_DEALLOCATE(diju_im)
         end if
         !Add FLL double-counting part
         if (ipert==0.and.pawu_new_algo.and.pawtab(itypat)%usepawu==5) then
           do ispden=1,nspdiag
             jrhoij=1
             do irhoij=1,pawrhoij(iatom)%nrhoijsel
               klmn=pawrhoij(iatom)%rhoijselect(irhoij)
               eldaumdc=eldaumdc+ro(1)*pawtab(itypat)%euij_fll(klmn)
               jrhoij=jrhoij+pawrhoij(iatom)%cplex
             end do
           end do
         end if
       end if
     end if
   end if

!  ========= Compute nuclear dipole moment energy contribution  ========
!  ==========================================================

!  Compute nuclear dipole contribution to Dij if necessary
   if (any(abs(nucdipmom(:,iatom))>tol8).and.ipert==0.and.ipositron/=1.and.option/=1) then
     if (paw_ij(iatom)%has_dijnd/=2) then
       call pawdijnd(cplex_dij,paw_ij(iatom)%dijnd,ndij,nucdipmom(:,iatom),pawrad(itypat),pawtab(itypat))
       paw_ij(iatom)%has_dijnd=2
     end if
   end if ! end compute dijnd if necessary

!  Compute contribution to on-site energy
   if (any(abs(nucdipmom(:,iatom))>tol8).and.ipert==0.and.ipositron/=1.and.option/=1) then
     ABI_CHECK(cplex==1,'BUG in pawdenpot: pawrhoij must be complex for ND moments!')
     jrhoij=2 !Select imaginary part of rhoij
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
       klmn=pawrhoij(iatom)%rhoijselect(irhoij)
!      Select imaginary part of dijnd (only nonzero part)
       kklmn=cplex_dij*(klmn-1)+2
!      Because dijnd involves no spin flips, following formula is correct for all values of nspden
       enucdip=enucdip-pawrhoij(iatom)%rhoijp(jrhoij,1)*paw_ij(iatom)%dijnd(kklmn,1)*pawtab(itypat)%dltij(klmn)
       jrhoij=jrhoij+cplex_rhoij
     end do
   end if

!  ========= Compute spin-orbit energy contribution  ========
!  ==========================================================

!  Compute spin-orbit contribution to Dij
   if (pawspnorb>0.and.ipert==0.and.ipositron/=1.and.(option/=2.or.cplex_rhoij==2)) then
     ABI_CHECK(cplex==1,'BUG in pawdenpot: Dij^SO not available for cplex=2!')
     call pawdijso(cplex,cplex_dij,paw_ij(iatom)%dijso,ndij,nspden,pawang,pawrad(itypat),pawtab(itypat), &
&                  pawxcdev,spnorbscl,paw_an(iatom)%vh1,paw_an(iatom)%vxc1)
     paw_ij(iatom)%has_dijso=2
     if (.not.keep_vhartree) then
       paw_an(iatom)%has_vhartree=0
       ABI_DEALLOCATE(paw_an(iatom)%vh1)
     end if
     if (temp_vxc) then
       paw_an(iatom)%has_vxc=0
       ABI_DEALLOCATE(paw_an(iatom)%vxc1)
     end if
   end if

!  Compute contribution to on-site energy
   if (option/=1.and.pawspnorb>0.and.ipert==0.and.ipositron/=1.and.cplex_rhoij==2) then
     ABI_CHECK(cplex==1,'BUG in pawdenpot: Epaw^SO not available for cplex=2!')
     ABI_CHECK(pawrhoij(iatom)%nspden==4,'BUG in pawdenpot: pawrhoij must have 4 components!')
     jrhoij=2 !Select imaginary part of rhoij
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
       klmn=pawrhoij(iatom)%rhoijselect(irhoij)
       kklmn=cplex_dij*(klmn-1)+1
       espnorb=espnorb-pawrhoij(iatom)%rhoijp(jrhoij,3)*paw_ij(iatom)%dijso(kklmn,3) &
&                     *pawtab(itypat)%dltij(klmn)
       if (cplex_dij==2) then
         kklmn=kklmn+1
         espnorb=espnorb-(pawrhoij(iatom)%rhoijp(jrhoij,2)*paw_ij(iatom)%dijso(kklmn,3) &
&              +half*pawrhoij(iatom)%rhoijp(jrhoij,4)*(paw_ij(iatom)%dijso(kklmn,1) &
&              -paw_ij(iatom)%dijso(kklmn,2))) &
&              *pawtab(itypat)%dltij(klmn)
       end if
       jrhoij=jrhoij+cplex_rhoij
     end do
   end if

!  === Compute 2nd part of local exact-exchange energy and potential  ===
!  ======================================================================

   if (pawtab(itypat)%useexexch>0.and.ipert==0.and.ipositron/=1) then

     ABI_CHECK(paw_ij(iatom)%nspden/=4,'BUG in pawdenpot: Local ex-exch. not implemented for nspden=4!')

     if (option<2) then
       call pawxpot(ndij,pawprtvol,pawrhoij(iatom),pawtab(itypat),paw_ij(iatom)%vpawx)
       paw_ij(iatom)%has_exexch_pot=2
     end if
     if (option/=1) then
       if (abs(pawprtvol)>=2) then
         write(msg, '(2a)' )ch10,'======= PAW local exact exchange terms (in Hartree) ===='
         call wrtout(std_out,  msg,'COLL')
         write(msg, '(2a,i4)' )ch10,' For Atom',iatom_tot
         call wrtout(std_out,  msg,'COLL')
       end if
       call pawxenergy(eexex,pawprtvol,pawrhoij(iatom),pawtab(itypat))
     end if

   end if ! useexexch

!  ==== Compute Fock Dij term and Fock energy terms ====
!  =====================================================

!  Computation of Fock contribution to Dij
   if (usefock==1) then

     ABI_CHECK(cplex==1,'BUG in pawdenpot: Dij^Fock not available for cplex=2!')

     if (ipositron/=1) then
!      Fock contribution to Dij
       ABI_ALLOCATE(dijfock_vv,(cplex_dij*lmn2_size,ndij))
       ABI_ALLOCATE(dijfock_cv,(cplex_dij*lmn2_size,ndij))
       call pawdijfock(cplex,cplex_dij,dijfock_vv,dijfock_cv,hyb_mixing_,hyb_mixing_sr_, &
&                      ndij,nspden,nsppol,pawrhoij(iatom),pawtab(itypat))
       paw_ij(iatom)%dijfock(:,:)=dijfock_vv(:,:)+dijfock_cv(:,:)

!      Fock contribution to energy

       if (option/=1) then
         if (cplex_dij==1) then
           do ispden=1,pawrhoij(iatom)%nspden
             jrhoij=1
             do irhoij=1,pawrhoij(iatom)%nrhoijsel
               klmn=pawrhoij(iatom)%rhoijselect(irhoij)
               ro(1)=pawrhoij(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
               efockdc=efockdc+ro(1)*half*dijfock_vv(klmn,ispden)
               efock=efock+ro(1)*(half*dijfock_vv(klmn,ispden)+dijfock_cv(klmn,ispden))
               jrhoij=jrhoij+cplex_rhoij
             end do
           end do
         else
           do ispden=1,pawrhoij(iatom)%nspden
             jrhoij=1
             do irhoij=1,pawrhoij(iatom)%nrhoijsel
               klmn=pawrhoij(iatom)%rhoijselect(irhoij);kklmn=2*klmn-1
               ro(1:2)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+1,ispden)*pawtab(itypat)%dltij(klmn)
               efockdc=efockdc+half*(ro(1)*dijfock_vv(kklmn,ispden)+ro(2)*dijfock_vv(kklmn+1,ispden))
               efock=efock+ro(1)*(half*dijfock_vv(kklmn,  ispden)+dijfock_cv(kklmn,  ispden)) &
&                         +ro(2)*(half*dijfock_vv(kklmn+1,ispden)+dijfock_cv(kklmn+1,ispden))
               jrhoij=jrhoij+cplex_rhoij
             end do
           end do
         end if
       end if

       ABI_DEALLOCATE(dijfock_vv)
       ABI_DEALLOCATE(dijfock_cv)

     else
!      Not Fock contribution for positron
       paw_ij(iatom)%dijfock(:,:)=zero
     end if
     paw_ij(iatom)%has_dijfock=2
   end if

!  ==========================================================
!  No more need of densities
   ABI_DEALLOCATE(rho1)
   ABI_DEALLOCATE(trho1)
   ABI_DEALLOCATE(nhat1)

!  === Compute the zero of the potentials if requested ==================
!  ======================================================================

   if (pawtab(itypat)%usepotzero==1 .AND. present(vpotzero)) then

     ! term 1 : beta
     vpotzero(1)=vpotzero(1)-pawtab(itypat)%beta/ucvol

     ! term 2 : \sum_ij rho_ij gamma_ij
     do ispden=1,nspdiag
       jrhoij=1
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         klmn=pawrhoij(iatom)%rhoijselect(irhoij)
         vpotzero(2) = vpotzero(2) - &
&         pawrhoij(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)*pawtab(itypat)%gammaij(klmn)/ucvol
         jrhoij=jrhoij+cplex_rhoij
       end do
     end do

   end if

!  =========== End loop on atoms ============================
!  ==========================================================

 end do

!========== Assemble "on-site" energy terms ===============
!==========================================================

 if (option/=1) then
   if (ipert==0) then
     epaw=e1xc+half*eh2+e1t10-exccore-etild1xc+eldaumdc+eexex+espnorb+efock+enucdip
     epawdc=e1xc-e1xcdc-half*eh2-exccore-etild1xc+etild1xcdc+eldaumdcdc-eexex-efockdc
   else
     epaw=e1xc-etild1xc+eh2+two*eldaumdc
     epawdc=zero
   end if
 end if

!========== Reduction in case of parallelism ==============
!==========================================================

 if (paral_atom) then
   if (option/=1)  then
     call timab(48,1,tsec)
     mpiarr=zero
     mpiarr(1)=compch_sph;mpiarr(2)=epaw;mpiarr(3)=epawdc
     if (ipositron/=0) then
       mpiarr(4)=electronpositron%e_paw
       mpiarr(5)=electronpositron%e_pawdc
     end if
     if (present(vpotzero)) then
       mpiarr(6)=vpotzero(1)
       mpiarr(7)=vpotzero(2)
     end if
     call xmpi_sum(mpiarr,my_comm_atom,ierr)
     compch_sph=mpiarr(1);epaw=mpiarr(2);epawdc=mpiarr(3)
     if (ipositron/=0) then
       electronpositron%e_paw=mpiarr(4)
       electronpositron%e_pawdc=mpiarr(5)
     end if
     if (present(vpotzero)) then
       vpotzero(1)=mpiarr(6)
       vpotzero(2)=mpiarr(7)
     end if
     call timab(48,2,tsec)
   end if
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 call timab(560,2,tsec)

 DBG_EXIT("COLL")

end subroutine pawdenpot
!!***

!----------------------------------------------------------------------

!!****f* m_paw_denpot/pawdensities
!! NAME
!! pawdensities
!!
!! FUNCTION
!! Compute PAW on-site densities (all-electron ,pseudo and compensation) for a given atom
!!
!! INPUTS
!!  cplex: if 1, on-site densities are REAL, if 2, COMPLEX (response function only)
!!  iatom=index of current atom (note: this is the absolute index, not the index on current proc)
!!  lm_size=number of (l,m) moments
!!  lmselectin(lm_size)=flags selecting the non-zero LM-moments of on-site densities
!!                      (value of these flags at input; must be .TRUE. for nzlmopt/=1)
!!  nspden=number of spin-density components
!!  nzlmopt=if -1, compute all LM-moments of densities (lmselectin=.true. forced)
!!                 initialize "lmselectout" (index of non-zero LM-moments of densities)
!!          if  0, compute all LM-moments of densities (lmselectin=.true. forced)
!!                 force "lmselectout" to .true. (index of non-zero LM-moments of densities)
!!          if  1, compute only non-zero LM-moments of densities (stored before in "lmselectin")
!!  one_over_rad2(mesh_size)= contains 1/r**2 for each point of the radial grid -optional argument-
!!  opt_compch=flag controlling the accumulation of compensation charge density moments
!!             inside PAW spheres (compch_sph)
!!  opt_dens=flag controlling which on-site density(ies) is (are) computed
!!           0: all on-site densities (all-electron, pseudo and compensation)
!!           1: all-electron and pseudo densities (no compensation)
!!           2: only all-electron density
!!  opt_l=controls which l-moment(s) contribute to the density:
!!        <0 : all l contribute
!!        >=0: only l=opt_l contributes
!!        Note: opt_l>=0 is only compatible with opt_dens=2
!!  opt_print=1 if the densities moments have to be printed out (only if pawprtvol>=2)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data (for the current atom type)
!!  pawrhoij <type(pawrhoij_type)>= paw rhoij occupancies and related data (for the current atom)
!!  pawtab <type(pawtab_type)>=paw tabulated starting data (for the current atom type)
!!
!! OUTPUT
!!  nhat1(cplex*mesh_size,lm_size,nspden)= compensation charge on-site density for current atom
!!  rho1(cplex*mesh_size,lm_size,nspden)= all electron on-site density for current atom
!!  trho1(cplex*mesh_size,lm_size,nspden)= pseudo on-site density for current atom
!!  ==== if nzlmopt/=1
!!    lmselectout(lm_size)=flags selecting the non-zero LM-moments of on-site densities
!!                         (value of these flags at output if updated, i.e. if nzlmopt<1)
!!
!!  SIDE EFFECTS
!!  ==== if opt_compch==1
!!    compch_sph=compensation charge integral inside spheres computed over spherical meshes
!!               updated with the contribution of current atom
!!
!! PARENTS
!!      make_efg_onsite,pawdenpot,pawdfptenergy,poslifetime,posratecore
!!
!! CHILDREN
!!      pawrad_deducer0,simp_gen,wrtout
!!
!! SOURCE

subroutine pawdensities(compch_sph,cplex,iatom,lmselectin,lmselectout,lm_size,nhat1,nspden,nzlmopt,&
&          opt_compch,opt_dens,opt_l,opt_print,pawang,pawprtvol,pawrad,pawrhoij,pawtab,rho1,trho1,&
&          one_over_rad2) ! optional


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdensities'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,iatom,lm_size,nspden,nzlmopt,opt_compch,opt_dens,opt_l,opt_print,pawprtvol
! jmb  real(dp),intent(out) :: compch_sph
 real(dp),intent(inout) :: compch_sph
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
 type(pawrhoij_type),intent(in) :: pawrhoij
 type(pawtab_type),intent(in) :: pawtab
!arrays
 logical,intent(in) :: lmselectin(lm_size)
 logical,intent(inout) :: lmselectout(lm_size)
 real(dp),intent(in),target,optional :: one_over_rad2(pawtab%mesh_size)
 real(dp),intent(out) :: nhat1(cplex*pawtab%mesh_size,lm_size,nspden*(1-((opt_dens+1)/2)))
 real(dp),intent(out) ::  rho1(cplex*pawtab%mesh_size,lm_size,nspden)
 real(dp),intent(out) :: trho1(cplex*pawtab%mesh_size,lm_size,nspden*(1-(opt_dens/2)))

!Local variables ---------------------------------------
!scalars
 integer :: dplex,ii,ilm,iplex,ir,irhoij,isel,ispden,jrhoij
 integer :: klm,klmn,kln,ll,lmax,lmin,mesh_size
 real(dp) :: m1,mt1,rdum
 character(len=500) :: msg
!arrays
 real(dp) :: compchspha(cplex),compchsphb(cplex),ro(cplex),ro_ql(cplex),ro_rg(cplex)
 real(dp),allocatable :: aa(:),bb(:)
 real(dp),pointer :: one_over_rad2_(:)

! *************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 if (opt_dens/=2.and.opt_l>=0) then
   msg='  opt_dens/=2 incompatible with opt_l>=0 !'
   MSG_BUG(msg)
 end if
 if(nzlmopt/=0.and.nzlmopt/=1.and.nzlmopt/=-1) then
   msg='  invalid value for variable "nzlmopt".'
   MSG_BUG(msg)
 end if
 if(nspden>pawrhoij%nspden) then
   msg='  nspden must be <= pawrhoij%nspden !'
   MSG_BUG(msg)
 end if
 if (cplex>pawrhoij%cplex) then
   msg='  cplex must be <= pawrhoij%cplex !'
   MSG_BUG(msg)
 end if
 if (nzlmopt/=1) then
   if (any(.not.lmselectin(1:lm_size))) then
     msg='  With nzlmopt/=1, lmselectin must be true !'
     MSG_BUG(msg)
   end if
 end if
 if (pawang%gnt_option==0) then
   msg='  pawang%gnt_option=0 !'
   MSG_BUG(msg)
 end if

!Various inits
 rho1=zero
 if (opt_dens <2) trho1=zero
 if (opt_dens==0) nhat1=zero
 mesh_size=pawtab%mesh_size;dplex=cplex-1
 if (nzlmopt<1) lmselectout(1:lm_size)=.true.
 if (present(one_over_rad2)) then
   one_over_rad2_ => one_over_rad2
 else
   ABI_ALLOCATE(one_over_rad2_,(mesh_size))
   one_over_rad2_(1)=zero
   one_over_rad2_(2:mesh_size)=one/pawrad%rad(2:mesh_size)**2
 end if

!===== Compute "on-site" densities (n1, ntild1, nhat1) =====
!==========================================================

 do ispden=1,nspden

!  -- Loop over ij channels (basis components)
   jrhoij=1
   do irhoij=1,pawrhoij%nrhoijsel
     klmn=pawrhoij%rhoijselect(irhoij)
     klm =pawtab%indklmn(1,klmn)
     kln =pawtab%indklmn(2,klmn)
     lmin=pawtab%indklmn(3,klmn)
     lmax=pawtab%indklmn(4,klmn)


!    Retrieve rhoij
     if (pawrhoij%nspden/=2) then
       ro(1:cplex)=pawrhoij%rhoijp(jrhoij:jrhoij+dplex,ispden)
     else
       if (ispden==1) then
         ro(1:cplex)=pawrhoij%rhoijp(jrhoij:jrhoij+dplex,1)&
&                   +pawrhoij%rhoijp(jrhoij:jrhoij+dplex,2)
       else if (ispden==2) then
         ro(1:cplex)=pawrhoij%rhoijp(jrhoij:jrhoij+dplex,1)
       end if
     end if
     ro(1:cplex)=pawtab%dltij(klmn)*ro(1:cplex)

!    First option: all on-site densities are computed (opt_dens==0)
!    --------------------------------------------------------------
     if (opt_dens==0) then
       do ll=lmin,lmax,2
         do ilm=ll**2+1,(ll+1)**2
           if (lmselectin(ilm)) then
             isel=pawang%gntselect(ilm,klm)
             if (isel>0) then
               ro_ql(1:cplex)=ro(1:cplex)*pawtab%qijl(ilm,klmn)
               ro_rg(1:cplex)=ro(1:cplex)*pawang%realgnt(isel)
!              == nhat1(r=0)
               nhat1(1:cplex,ilm,ispden)=nhat1(1:cplex,ilm,ispden) &
&               +ro_ql(1:cplex)*pawtab%shapefunc(1,ll+1)
!              == rho1(r>0), trho1(r>0), nhat1(r>0)
               do ir=2,mesh_size
                 rho1(cplex*ir-dplex:ir*cplex,ilm,ispden) =rho1(cplex*ir-dplex:ir*cplex,ilm,ispden)&
&                 +ro_rg(1:cplex)*pawtab%phiphj  (ir,kln)*one_over_rad2_(ir)
                 trho1(cplex*ir-dplex:ir*cplex,ilm,ispden)=trho1(cplex*ir-dplex:ir*cplex,ilm,ispden)&
&                 +ro_rg(1:cplex)*pawtab%tphitphj(ir,kln)*one_over_rad2_(ir)
                 nhat1(cplex*ir-dplex:ir*cplex,ilm,ispden)=nhat1(cplex*ir-dplex:ir*cplex,ilm,ispden)&
&                 +ro_ql(1:cplex)*pawtab%shapefunc(ir,ll+1)
               end do
             end if
           end if
         end do  ! End loops over ll,lm
       end do

!      2nd option: AE and pseudo densities are computed (opt_dens==1)
!      --------------------------------------------------------------
     else if (opt_dens==1) then
       do ll=lmin,lmax,2
         do ilm=ll**2+1,(ll+1)**2
           if (lmselectin(ilm)) then
             isel=pawang%gntselect(ilm,klm)
             if (isel>0) then
               ro_rg(1:cplex)=ro(1:cplex)*pawang%realgnt(isel)
!              == rho1(r>0), trho1(r>0)
               do ir=2,mesh_size
                 rho1(cplex*ir-dplex:ir*cplex,ilm,ispden) =rho1(cplex*ir-dplex:ir*cplex,ilm,ispden)&
&                 +ro_rg(1:cplex)*pawtab%phiphj  (ir,kln)*one_over_rad2_(ir)
                 trho1(cplex*ir-dplex:ir*cplex,ilm,ispden)=trho1(cplex*ir-dplex:ir*cplex,ilm,ispden)&
&                 +ro_rg(1:cplex)*pawtab%tphitphj(ir,kln)*one_over_rad2_(ir)
               end do
             end if
           end if
         end do  ! End loops over ll,lm
       end do

!      3rd option: only all-electron on-site density is computed (opt_dens==2)
!      -----------------------------------------------------------------------
     else if (opt_dens==2) then
       if (opt_l<0.or.(pawtab%indklmn(3,klmn)==0.and.pawtab%indklmn(4,klmn)==2*opt_l)) then
         do ll=lmin,lmax,2
           do ilm=ll**2+1,(ll+1)**2
             if (lmselectin(ilm)) then
               isel=pawang%gntselect(ilm,klm)
               if (isel>0) then
                 ro_rg(1:cplex)=ro(1:cplex)*pawang%realgnt(isel)
!                == rho1(r>0)
                 do ir=2,mesh_size
                   rho1(cplex*ir-dplex:ir*cplex,ilm,ispden) =rho1(cplex*ir-dplex:ir*cplex,ilm,ispden)&
&                   +ro_rg(1:cplex)*pawtab%phiphj(ir,kln)*one_over_rad2_(ir)
                 end do
               end if
             end if
           end do  ! End loops over ll, lm
         end do
       end if
     end if

!    -- End loop over ij channels
     jrhoij=jrhoij+pawrhoij%cplex
   end do

!  Scale densities with 1/r**2 and compute rho1(r=0) and trho1(r=0)
   if (cplex==2)  then
     ABI_ALLOCATE(aa,(5))
     ABI_ALLOCATE(bb,(5))
   end if
   if (opt_dens==0.or.opt_dens==1) then
     do ll=0,pawtab%lcut_size-1
       do ilm=ll**2+1,(ll+1)**2
         if (lmselectin(ilm)) then
           if (cplex==1) then
             call pawrad_deducer0(rho1 (:,ilm,ispden),mesh_size,pawrad)
             call pawrad_deducer0(trho1(:,ilm,ispden),mesh_size,pawrad)
           else
             do ii=0,1
               do ir=2,5
                 aa(ir)=rho1 (2*ir-ii,ilm,ispden)
                 bb(ir)=trho1(2*ir-ii,ilm,ispden)
               end do
               call pawrad_deducer0(aa,5,pawrad)
               call pawrad_deducer0(bb,5,pawrad)
               rho1 (2-ii,ilm,ispden)=aa(1)
               trho1(2-ii,ilm,ispden)=bb(1)
             end do
           end if
         end if
       end do
     end do
   else
     do ll=0,pawtab%lcut_size-1
       do ilm=ll**2+1,(ll+1)**2
         if (lmselectin(ilm)) then
           if (cplex==1) then
             call pawrad_deducer0(rho1(:,ilm,ispden),mesh_size,pawrad)
           else
             do ii=0,1
               do ir=2,5
                 aa(ir)=rho1 (2*ir-ii,ilm,ispden)
               end do
               call pawrad_deducer0(aa,5,pawrad)
               rho1(2-ii,ilm,ispden)=aa(1)
             end do
           end if
         end if
       end do
     end do
   end if
   if (cplex==2)  then
     ABI_DEALLOCATE(aa)
     ABI_DEALLOCATE(bb)
   end if

!  -- Test moments of densities and store non-zero ones
   if (nzlmopt==-1) then
     do ll=0,pawtab%lcut_size-1
       do ilm=ll**2+1,(ll+1)**2
         m1=zero;mt1=zero
         if (cplex==1) then
           m1=maxval(abs(rho1 (1:mesh_size,ilm,ispden)))
           if (opt_dens<2) mt1=maxval(abs(trho1(1:mesh_size,ilm,ispden)))
         else
           do ir=1,mesh_size
             rdum=sqrt(rho1(2*ir-1,ilm,ispden)**2+rho1(2*ir,ilm,ispden)**2)
             m1=max(m1,rdum)
           end do
           if (opt_dens<2) then
             do ir=1,mesh_size
               rdum=sqrt(trho1(2*ir-1,ilm,ispden)**2+trho1(2*ir,ilm,ispden)**2)
               mt1=max(mt1,rdum)
             end do
           end if
         end if
         if (ispden==1) then
           if ((ilm>1).and.(m1<tol16).and.(mt1<tol16)) then
             lmselectout(ilm)=.false.
           end if
         else if (.not.(lmselectout(ilm))) then
           lmselectout(ilm)=((m1>=tol16).or.(mt1>=tol16))
         end if
       end do
     end do
   end if

!  -- Compute integral of (n1-tn1) inside spheres
   if (opt_compch==1.and.ispden==1.and.opt_dens<2) then
     ABI_ALLOCATE(aa,(mesh_size))
     aa(1:mesh_size)=(rho1(1:mesh_size,1,1)-trho1(1:mesh_size,1,1)) &
&     *pawrad%rad(1:mesh_size)**2
! jmb
!    compchspha(1)=zero
     call simp_gen(compchspha(1),aa,pawrad)
     compch_sph=compch_sph+compchspha(1)*sqrt(four_pi)
     ABI_DEALLOCATE(aa)
   end if

!  -- Print out moments of densities (if requested)
   if (abs(pawprtvol)>=2.and.opt_print==1.and.opt_dens<2) then
     ABI_ALLOCATE(aa,(cplex*mesh_size))
     ABI_ALLOCATE(bb,(cplex*mesh_size))
     if (opt_dens==0) then
       write(msg,'(2a,i3,a,i1,3a)') ch10, &
&       ' Atom ',iatom,' (ispden=',ispden,'):',ch10,&
&       '  ******* Moment of (n1-tn1)  ** Moment of (n1-tn1-nhat1)'
     else
       write(msg,'(2a,i3,a,i1,3a)') ch10, &
&       ' Atom ',iatom,' (ispden=',ispden,'):',ch10,&
&       '  ******* Moment of (n1-tn1)'
     end if
     call wrtout(std_out,msg,'PERS')
     do ll=0,pawtab%lcut_size-1
       do ilm=ll**2+1,(ll+1)**2
         if (lmselectin(ilm)) then
           do iplex=1,cplex
             if (opt_dens==0) then
               do ir=1,mesh_size
                 ii=cplex*(ir-1)+iplex
                 ro(1)=pawrad%rad(ir)**(2+ll)
                 aa(ir)=ro(1)*(rho1(ii,ilm,ispden)-trho1(ii,ilm,ispden))
                 bb(ir)=ro(1)*nhat1(ii,ilm,ispden)
               end do
               call simp_gen(compchspha(iplex),aa,pawrad)
               call simp_gen(compchsphb(iplex),bb,pawrad)
             else
               do ir=1,mesh_size
                 ii=cplex*(ir-1)+iplex
                 ro(1)=pawrad%rad(ir)**(2+ll)
                 aa(ir)=ro(1)*(rho1(ii,ilm,ispden)-trho1(ii,ilm,ispden))
               end do
               call simp_gen(compchspha(iplex),aa,pawrad)
             end if
           end do
           if (opt_dens==0) then
             if (cplex==1) then
               write(msg,'(3x,a,2i2,2(a,es14.7))') &
&               'l,m=',ll,ilm-(ll**2+ll+1),': M=',compchspha(1),&
&               ' **    M=',compchspha(1)-compchsphb(1)
             else
               write(msg,'(3x,a,2i2,2(a,2es14.7))') &
&               'l,m=',ll,ilm-(ll**2+ll+1),': M=',compchspha(1:2),&
&               ' **    M=',compchspha(1:2)-compchsphb(1:2)
             end if
           else
             if (cplex==1) then
               write(msg,'(3x,a,2i2,a,es14.7)') &
&               'l,m=',ll,ilm-(ll**2+ll+1),': M=',compchspha(1)
             else
               write(msg,'(3x,a,2i2,a,2es14.7)') &
&               'l,m=',ll,ilm-(ll**2+ll+1),': M=',compchspha(1:2)
             end if
           end if
           call wrtout(std_out,msg,'PERS')
         end if
       end do
     end do
     ABI_DEALLOCATE(aa)
     ABI_DEALLOCATE(bb)
   end if

!  ----- End loop over spin components
 end do

 if (.not.present(one_over_rad2))  then
   ABI_DEALLOCATE(one_over_rad2_)
 end if

 DBG_EXIT("COLL")

end subroutine pawdensities
!!***

!----------------------------------------------------------------------

!!****f* m_paw_denpot/paw_mknewh0
!! NAME
!! paw_mknewh0
!!
!! FUNCTION
!! Calculates the new bare PAW Hamiltonian in the case of quasi-particle self-consistent GW calculations.
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nspden=number of spin-density components
!!  nfftf=(effective) number of FFT grid points (for this proc) for the "fine" grid
!!  pawspnorb=flag: 1 if spin-orbit coupling is activated
!!  pawprtvol=control print volume and debugging output for PAW
!!  Cryst<crystal_t>=Info on unit cell and its symmetries
!!  Pawtab(ntypat*usepaw)<type(pawtab_type)>=paw tabulated starting data
!!  Paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  Pawang<type(pawang_type)>=paw angular mesh and related data
!!  Pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  vxc(nfftf,nspden)=exchange-correlation potential
!!  vxc_val(nfftf,nspden)=valence only exchange-correlation potential
!!  vtrial(nfftf,nspden)=potential (Hartree+XC+loc)
!!
!! SIDE EFFECTS
!!  Paw_ij(natom*usepaw)<Paw_ij_type)>=paw arrays given on (i,j) channels
!!     At output: new value for Paw_ij()%dij
!!
!! PARENTS
!!      calc_vhxc_me
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,pawgylm,symdij,symdij_all,wrtout
!!
!! SOURCE

subroutine paw_mknewh0(my_natom,nsppol,nspden,nfftf,pawspnorb,pawprtvol,Cryst,&
&          Pawtab,Paw_an,Paw_ij,Pawang,Pawfgrtab,vxc,vxc_val,vtrial,&
&          mpi_atmtab,comm_atom) ! optional arguments (parallelism)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_mknewh0'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom,nsppol,nspden,nfftf,pawprtvol,pawspnorb
 integer,optional,intent(in) :: comm_atom
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: vxc(nfftf,nspden),vxc_val(nfftf,nspden),vtrial(nfftf,nspden)
 type(crystal_t),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(Pawtab_type),target,intent(in) :: Pawtab(Cryst%ntypat)
 type(Paw_an_type),intent(in) :: Paw_an(my_natom)
 type(Paw_ij_type),intent(inout) :: Paw_ij(my_natom)
 type(Pawfgrtab_type),intent(inout) :: Pawfgrtab(my_natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: ipert0=0
 integer :: iat,iat_tot,idij,cplex_rf,ndij,option_dij
 integer :: itypat,lmn_size,j0lmn,jlmn,ilmn,klmn,klmn1,klm
 integer :: lmin,lmax,mm,isel,lm_size,lmn2_size,my_comm_atom,cplex_dij
 integer :: ils,ilslm,ic,lm0
 integer :: nsploop,is2fft
 real(dp) :: gylm,qijl
 logical :: ltest,my_atmtab_allocated,paral_atom
 character(len=500) :: msg
!arrays
 integer, pointer :: indklmn_(:,:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: rdum(1),rdum2(1)
 real(dp),allocatable :: prod_hloc(:,:),prodhxc_core(:,:)
 real(dp),allocatable :: dijhl_hat(:,:),dijhmxc_val(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call wrtout(std_out,'Assembling PAW strengths for the bare Hamiltonian','COLL')

!== Set up parallelism over atoms ===
 paral_atom=(present(comm_atom).and.(my_natom/=Cryst%natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,Cryst%natom,my_natom_ref=my_natom)

 if (my_natom>0) then

!  === Test if required pointers in paw_ij are allocated ===
   ltest = (allocated(Paw_ij(1)%dijxc).and.allocated(Paw_ij(1)%dijxc_val) )
   ABI_CHECK(ltest,'dijxc or dijxc_val not calculated')

   ltest=(allocated(Paw_ij(1)%dijhat)) !.and.Paw_ij(1)%has_dijhat==2)
   ABI_CHECK(ltest,'dijhat not calculated')

   ltest=(allocated(Paw_ij(1)%dijhartree)) !.and.Paw_ij(1)%has_dijhartree==2)
   ABI_CHECK(ltest,'dijhartree not calculated')

   if (ANY(Pawtab(:)%usepawu>0)) then
     do iat=1,my_natom
       iat_tot=iat;if (paral_atom) iat_tot=my_atmtab(iat)
       itypat=Cryst%typat(iat_tot)
       if (Pawtab(itypat)%usepawu>0) then
         ltest=(allocated(Paw_ij(iat)%dijU) ) !.and.Paw_ij(iat)%has_dijU==2)
         write(msg,'(a,i3,a)')" For atom no. ",iat," %dijU(iat) has not been calculated."
         ABI_CHECK(ltest,msg)
       end if
     end do
   end if

   if (pawspnorb>0) then
     do iat=1,my_natom
       ltest=(allocated(Paw_ij(iat)%dijso) ) !.and.Paw_ij(iat)%has_dijso==2)
       write(msg,'(a,i3,a)')" For atom no. ",iat," %dijso(iat) has not been calculated."
       ABI_CHECK(ltest,msg)
     end do
   end if
 end if ! my_natom>0

!== Construct the new PAW H0 Hamiltonian ===
 do iat=1,my_natom
   iat_tot=iat;if (paral_atom) iat_tot=my_atmtab(iat)

   itypat    = Cryst%typat(iat_tot)
   lmn_size  = Pawtab(itypat)%lmn_size
   lmn2_size = Pawtab(itypat)%lmn2_size
   lm_size   = Paw_an(iat)%lm_size
   cplex_rf  = Paw_ij(iat)%cplex_rf
   cplex_dij = Paw_ij(iat)%cplex_dij
   ndij      = Paw_ij(iat)%ndij

   ABI_CHECK(cplex_rf==1,'cplex_rf/=1 not implemented')
   ABI_CHECK(cplex_dij==1,'cplex_dij/=1 not implemented')
!
!  Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
   if (Pawfgrtab(iat)%gylm_allocated==0) then
     if (allocated(Pawfgrtab(iat)%gylm))  then
       ABI_DEALLOCATE(Pawfgrtab(iat)%gylm)
     end if
     ABI_ALLOCATE(Pawfgrtab(iat)%gylm,(Pawfgrtab(iat)%nfgd,lm_size))
     Pawfgrtab(iat)%gylm_allocated=2

     call pawgylm(Pawfgrtab(iat)%gylm,rdum,rdum2,lm_size,&
&     Pawfgrtab(iat)%nfgd,1,0,0,Pawtab(itypat),Pawfgrtab(iat)%rfgd)
   end if

!  === Calculate LM contribution to dijhmxc_val for this atom ===
!  * Dijxc contains also the Hat term on the FFT mesh while Dijxc_val does not
!  contain neither the hat term nor the LM sum of onsite terms (they should cancel each other)
!  FIXME change paw_dij,  otherwise I miss tnc in vxc
!  * prodhxc_core is used to assemble $\int g_l Ylm (vtrial - vxc_val[tn+nhat] dr$ on the FFT mesh ===
!  * The following quantities do not depend on ij
   ABI_ALLOCATE(prod_hloc   ,(lm_size,ndij))
   ABI_ALLOCATE(prodhxc_core,(lm_size,ndij))
   prod_hloc   =zero
   prodhxc_core=zero
   do idij=1,ndij
     do ilslm=1,lm_size
       do ic=1,Pawfgrtab(iat)%nfgd
         is2fft=Pawfgrtab(iat)%ifftsph(ic)
         gylm=Pawfgrtab(iat)%gylm(ic,ilslm)
         prod_hloc (ilslm,idij)=prod_hloc (ilslm,idij) + (vtrial(is2fft,idij)-vxc(is2fft,idij))*gylm
!        prodhxc_core(ilslm,idij)=prodhxc_core(ilslm,idij) + (vxc_val(is2fft,idij))*gylm
         prodhxc_core(ilslm,idij)=prodhxc_core(ilslm,idij) + (vtrial(is2fft,idij)-vxc_val(is2fft,idij))*gylm
       end do
     end do
   end do !idij

!  === Assembly the "Hat" contribution for this atom ====
   ABI_ALLOCATE(dijhl_hat  ,(cplex_dij*lmn2_size,ndij))
   ABI_ALLOCATE(dijhmxc_val,(cplex_dij*lmn2_size,ndij))
   dijhl_hat  =zero
   dijhmxc_val=zero
   indklmn_ => Pawtab(itypat)%indklmn(1:6,1:lmn2_size)

   do idij=1,ndij
     do klmn=1,lmn2_size
       klm =indklmn_(1,klmn)
       lmin=indklmn_(3,klmn)
       lmax=indklmn_(4,klmn)

!      === $\sum_lm q_ij^l prod* for each idij$ ===
       do ils=lmin,lmax,2
         lm0=ils**2+ils+1
         do mm=-ils,ils
           ilslm=lm0+mm
           isel=Pawang%gntselect(lm0+mm,klm)
           if (isel>0) then
             qijl=Pawtab(itypat)%qijl(ilslm,klmn)
             dijhl_hat  (klmn,idij)=dijhl_hat  (klmn,idij) +  prod_hloc (ilslm,idij)*qijl
             dijhmxc_val(klmn,idij)=dijhmxc_val(klmn,idij) +prodhxc_core(ilslm,idij)*qijl
           end if
         end do
       end do
     end do
   end do

   ABI_DEALLOCATE(prod_hloc)
   ABI_DEALLOCATE(prodhxc_core)

!  * Normalization factor due to integration on the FFT mesh
   dijhl_hat  = dijhl_hat  *Cryst%ucvol/DBLE(nfftf)
   dijhmxc_val= dijhmxc_val*Cryst%ucvol/DBLE(nfftf)

!  === Now assembly the bare Hamiltonian ===
!  * Loop over density components overwriting %dij
   nsploop=nsppol; if (Paw_ij(iat)%ndij==4) nsploop=4

   do idij=1,nsploop
     klmn1=1

     do jlmn=1,lmn_size
       j0lmn=jlmn*(jlmn-1)/2
       do ilmn=1,jlmn
         klmn=j0lmn+ilmn

!        The following gives back the input dij.
!        since dijxc contains the hat term done on the FFT mesh
         if (.FALSE.) then
           Paw_ij(iat)%dij(klmn,idij) =        &
&           Pawtab(itypat)%dij0    (klmn)      &
&           +Paw_ij(iat)%dijhartree(klmn)      &
&           +Paw_ij(iat)%dijxc     (klmn,idij) &
&           +dijhl_hat   (klmn,idij)

         else
!          === Make nonlocal part of h0 removing the valence contribution ===
!          Remeber that XC contains already the Hat contribution
           Paw_ij(iat)%dij(klmn,idij) =        &
&           Pawtab(itypat)%dij0      (klmn)    &
&           +Paw_ij(iat)%dijhartree(klmn)      &
&           +Paw_ij(iat)%dijxc     (klmn,idij) &  ! 2 lines to get the d1-dt1 XC core contribution + XC hat (core+val)
&          -Paw_ij(iat)%dijxc_val (klmn,idij) &  ! I suppose that the "hat" term on the FFT mesh in included in both.
&          +dijhmxc_val(klmn,idij)               ! Local + Hartree - XC val contribution to the "hat" term.

!          Add the U contribution to the
!          if (.FALSE. .and. Pawtab(itypat)%usepawu>0) then
           if (.TRUE. .and. Pawtab(itypat)%usepawu>0) then
             Paw_ij(iat)%dij(klmn,idij) = Paw_ij(iat)%dij(klmn,idij) + Paw_ij(iat)%dijU(klmn,idij)
           end if
         end if
!        TODO dijso, dijU, vpawx?
!        Just to be consistent, update some values.
!$Paw_ij(iat)%dijhat(klmn,idij)=Paw_ij(iat)%dijhat(klmn,idij)-dijhmxc_val(klmn,idij)

       end do !ilmn
     end do !jlmn
   end do !idij

!  this is to be consistent?
!  deallocate(Paw_ij(iat)%dijvxc_val)
   ABI_DEALLOCATE(dijhl_hat)
   ABI_DEALLOCATE(dijhmxc_val)
 end do !iat

!=== Symmetrize total Dij ===
 option_dij=0 ! For total Dij.
#if 0
 if (paral_atom) then
   call symdij(Cryst%gprimd,Cryst%indsym,ipert0,my_natom,Cryst%natom,Cryst%nsym,Cryst%ntypat,option_dij,&
&   Paw_ij,Pawang,pawprtvol,Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec,&
&   comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
 else
   call symdij(Cryst%gprimd,,Cryst%indsym,ipert0,my_natom,Cryst%natom,Cryst%nsym,Cryst%ntypat,option_dij,&
&   Paw_ij,Pawang,pawprtvol,Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec)
 end if
#else
 if (paral_atom) then
   call symdij_all(Cryst%gprimd,Cryst%indsym,ipert0,my_natom,Cryst%natom,Cryst%nsym,Cryst%ntypat,&
&   Paw_ij,Pawang,pawprtvol,Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec,&
&   comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
 else
   call symdij_all(Cryst%gprimd,Cryst%indsym,ipert0,my_natom,Cryst%natom,Cryst%nsym,Cryst%ntypat,&
&   Paw_ij,Pawang,pawprtvol,Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec)
 end if
#endif

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine paw_mknewh0
!!***

!----------------------------------------------------------------------

END MODULE m_paw_denpot
!!***

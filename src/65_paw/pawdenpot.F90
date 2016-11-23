!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawdenpot
!! NAME
!! pawdenpot
!!
!! FUNCTION
!! Compute different (PAW) energies, densities and potentials (or potential-like quantities)
!! inside PAW spheres
!! Can also compute first-order densities potentials and second-order energies (RF calculations).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!!    paw_ij(my_natom)%dijso(cplex_dij*lmn2_size,nspden)=spin-orbit contribution to dij
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
!!      pawdijnd,pawdijso,pawrad_deducer0,pawuenergy,pawxc,pawxc3,pawxcm
!!      pawxcm3,pawxcmpositron,pawxcpositron,pawxenergy,pawxpot,poisson
!!      setnoccmmp,timab,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawdenpot(compch_sph,epaw,epawdc,ipert,ixc,&
& my_natom,natom,nspden,ntypat,nucdipmom,nzlmopt,option,paw_an,paw_an0,&
& paw_ij,pawang,pawprtvol,pawrad,pawrhoij,pawspnorb,pawtab,pawxcdev,spnorbscl,xclevel,xc_denpos,ucvol,znucl,&
& electronpositron,mpi_atmtab,comm_atom,vpotzero,hyb_mixing,hyb_mixing_sr) ! optional arguments

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_pawang,  only: pawang_type
 use m_pawrad,  only: pawrad_type, pawrad_deducer0, poisson, simp_gen
 use m_pawtab,  only: pawtab_type
 use m_paw_an, only : paw_an_type
 use m_paw_ij, only : paw_ij_type
 use m_pawrhoij,only: pawrhoij_type
 use m_pawdij,  only: pawdijhartree, pawdijnd, pawdijso, pawxpot ,pawdijfock
 use m_pawxc,   only: pawxc, pawxc3, pawxcm, pawxcm3, pawxcpositron, pawxcmpositron
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab
 use m_electronpositron, only: electronpositron_type,electronpositron_calctype

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdenpot'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_65_paw, except_this_one => pawdenpot
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
 integer :: cplex,cplex_dij,has_kxc,iatom,iatom_tot,idum,ierr,ipositron,irhoij,ispden,itypat,itypat0
 integer :: jrhoij,kklmn,klmn,lm_size,lmn2_size,mesh_size,my_comm_atom,ndij,nkxc1,nspdiag,nsppol,opt_compch
 integer :: usecore,usepawu,usetcore,usexcnhat,usefock
 logical :: keep_vhartree,my_atmtab_allocated,need_kxc,paral_atom,temp_vxc
 real(dp) :: e1t10,e1xc,e1xcdc,efock,efockdc,eexc,eexcdc,eexdctemp
 real(dp) :: eexc_val,eexcdc_val,eexex,eexexdc,eextemp,eh2
 real(dp) :: eldaumdc,eldaumdcdc,enucdip,espnorb,etild1xc,etild1xcdc
 real(dp) :: exccore,exchmix,hyb_mixing_,hyb_mixing_sr_,rdum
 character(len=3) :: pertstrg
 character(len=500) :: msg
!arrays
 integer :: idum1(0),idum3(0,0,0)
 integer,pointer :: my_atmtab(:)
 logical,allocatable :: lmselect_cur(:),lmselect_cur_ep(:),lmselect_ep(:),lmselect_tmp(:)
 real(dp) :: ro(2),mpiarr(7),tsec(2)
 real(dp),allocatable :: dij_ep(:),dijfock_vv(:,:),dijfock_cv(:,:),one_over_rad2(:),kxc_tmp(:,:,:),nhat1(:,:,:),nhat1_ep(:,:,:)
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

 if (my_natom>0) then
   if(paw_ij(1)%has_dijhartree==0.and. .not.(ipert==natom+1.or.ipert==natom+10.or.ipert==natom+11)) then
     msg='dijhartree must be allocated !'
     MSG_BUG(msg)
   end if
   if (paw_ij(1)%cplex/=paw_an(1)%cplex) then
     msg='paw_ij()%cplex and paw_an()%cplex must be equal !'
     MSG_BUG(msg)
   end if
   if (pawrhoij(1)%cplex<paw_an(1)%cplex) then
     msg='pawrhoij()%cplex must be >=paw_an()%cplex  !'
     MSG_BUG(msg)
   end if
   if (ipert<=0.and.paw_ij(1)%cplex/=1) then
     msg='cplex must be 1 for GS calculations !'
     MSG_BUG(msg)
   end if
   if (ipert>0.and.(ipert<=natom.or.ipert==natom+2).and.paw_an0(1)%has_kxc/=2) then
     msg='XC kernels for ground state must be in memory !'
     MSG_BUG(msg)
   end if
   if(paw_an(1)%has_vxc==0.and.(option==0.or.option==1).and. &
&   .not.(ipert==natom+1.or.ipert==natom+10.or.ipert==natom+11)) then
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
 if (ipert==natom+1.or.ipert==natom+10.or.ipert==natom+11) then
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
 usepawu=maxval(pawtab(1:ntypat)%usepawu)
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
   has_kxc=paw_an(iatom)%has_kxc;need_kxc=(has_kxc==1)
   cplex=1;if (ipert>0) cplex=pawrhoij(iatom)%cplex
   cplex_dij=paw_ij(iatom)%cplex_dij
   ndij=paw_ij(iatom)%ndij

!  Allocations of "on-site" densities
   ABI_ALLOCATE(rho1 ,(cplex*mesh_size,lm_size,nspden))
   ABI_ALLOCATE(trho1,(cplex*mesh_size,lm_size,nspden))
   ABI_ALLOCATE(nhat1,(cplex*mesh_size,lm_size,nspden*usexcnhat))
   rho1(:,:,:)=zero;trho1(:,:,:)=zero;nhat1(:,:,:)=zero
   if (ipositron/=0) then ! Additional allocation for the electron-positron case
     ABI_ALLOCATE(rho1_ep ,(cplex*mesh_size,lm_size,nspden))
     ABI_ALLOCATE(trho1_ep,(cplex*mesh_size,lm_size,nspden))
     ABI_ALLOCATE(nhat1_ep,(cplex*mesh_size,lm_size,nspden*usexcnhat))
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
&   pawrhoij(iatom)%cplex==2.and.paw_an(iatom)%has_vxc==0) then
!    these should already be allocated in paw_an_init!
     if (pawxcdev==0)then
       if (allocated(paw_an(iatom)%vxc1))  then
         ABI_DEALLOCATE(paw_an(iatom)%vxc1)
       end if
       ABI_ALLOCATE(paw_an(iatom)%vxc1,(cplex*mesh_size,paw_an(iatom)%angl_size,nspden))
     end if
     if (pawxcdev/=0 ) then
       if (allocated(paw_an(iatom)%vxc1))  then
         ABI_DEALLOCATE(paw_an(iatom)%vxc1)
       end if
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
&   nhat1,nspden,nzlmopt,opt_compch,1-usexcnhat,-1,1,pawang,pawprtvol,pawrad(itypat),&
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
&     lm_size,nhat1_ep,nspden,nzlmopt,0,1-usexcnhat,-1,0,pawang,0,pawrad(itypat),&
&     electronpositron%pawrhoij_ep(iatom),pawtab(itypat),&
&     rho1_ep,trho1_ep,one_over_rad2=one_over_rad2)

     if (nzlmopt<1) electronpositron%lmselect_ep(1:lm_size,iatom)=lmselect_ep(1:lm_size)
     ABI_DEALLOCATE(lmselect_ep)
     ABI_DEALLOCATE(lmselect_cur_ep)
   end if

!  =========== Compute XC potentials and energies ===========
!  ==========================================================

!  Temporary storage
   nkxc1=0;if (paw_an(iatom)%has_kxc/=0) nkxc1=paw_an(iatom)%nkxc1
   if (pawxcdev/=0) then
     ABI_ALLOCATE(vxc_tmp,(cplex*mesh_size,lm_size,nspden))
     if (need_kxc)  then
       ABI_ALLOCATE(kxc_tmp,(mesh_size,lm_size,nkxc1))
     end if
   end if
   if (pawxcdev==0) then
     ABI_ALLOCATE(vxc_tmp,(cplex*mesh_size,pawang%angl_size,nspden))
     vxc_tmp(:,:,:)=zero
     if (need_kxc)  then
       ABI_ALLOCATE(kxc_tmp,(mesh_size,pawang%angl_size,nkxc1))
     end if
   end if
   idum=0
   if (.not.allocated(kxc_tmp))  then
     ABI_ALLOCATE(kxc_tmp,(0,0,0))
   end if

!  ===== Vxc1 term =====
   if (ipositron/=1) then
     if (pawxcdev/=0) then
       if (ipert==0) then
         call pawxcm(pawtab(itypat)%coredens,eexc,eexcdc,idum,ixc,kxc_tmp,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,nkxc1,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),pawxcdev,rho1,usecore,0,vxc_tmp,xclevel,xc_denpos)
       else
         call pawxcm3(pawtab(itypat)%coredens,cplex,cplex,eexc,ixc,paw_an0(iatom)%kxc1,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,paw_an0(iatom)%nkxc1,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),rho1,usecore,0,vxc_tmp,xclevel)
         eexcdc=zero
       end if
     else
       if (ipert==0) then
         call pawxc(pawtab(itypat)%coredens,eexc,eexcdc,ixc,kxc_tmp,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,nkxc1,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),rho1,usecore,0,vxc_tmp,xclevel,xc_denpos)
       else
         call pawxc3(pawtab(itypat)%coredens,cplex,cplex,eexc,ixc,paw_an0(iatom)%kxc1,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,paw_an0(iatom)%nkxc1,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),rho1,usecore,0,vxc_tmp,xclevel)
         eexcdc=zero
       end if
     end if

     if (option/=1) then
       e1xc=e1xc+eexc
       e1xcdc=e1xcdc+eexcdc
     end if
     if (option<2.or.temp_vxc) paw_an(iatom)%vxc1(:,:,:)=vxc_tmp(:,:,:)
     if (need_kxc.and.nkxc1>0) paw_an(iatom)%kxc1(:,:,:)=kxc_tmp(:,:,:)
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
           !ABI_ALLOCATE(kxc_tmp,(0,0,0))
         call pawxcm(pawtab(itypat)%tcoredens(:,1),&
&         eexc,eexcdc,idum,ixc,kxc_tmp,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,nkxc1,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),pawxcdev,trho1,usetcore,2*usexcnhat,vxc_tmp,xclevel,xc_denpos)
           !ABI_DEALLOCATE(kxc_tmp)
       else
         call pawxcm3(pawtab(itypat)%tcoredens(:,1),&
&         cplex,cplex,eexc,ixc,paw_an0(iatom)%kxct1,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,paw_an0(iatom)%nkxc1,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),trho1,usetcore,2*usexcnhat,vxc_tmp,xclevel)
         eexcdc=zero
       end if
     else
       if (ipert==0) then
         call pawxc(pawtab(itypat)%tcoredens(:,1),&
&         eexc,eexcdc,ixc,kxc_tmp,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,nkxc1,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),trho1,usetcore,2*usexcnhat,vxc_tmp,xclevel,xc_denpos)
       else
         call pawxc3(pawtab(itypat)%tcoredens(:,1),&
&         cplex,cplex,eexc,ixc,paw_an0(iatom)%kxct1,lm_size,&
&         paw_an(iatom)%lmselect,nhat1,paw_an0(iatom)%nkxc1,mesh_size,nspden,option,&
&         pawang,pawrad(itypat),trho1,usetcore,2*usexcnhat,vxc_tmp,xclevel)
         eexcdc=zero
       end if
     end if

     if (option/=1) then
       etild1xc=etild1xc+eexc
       etild1xcdc=etild1xcdc+eexcdc
     end if
     if (option<2) paw_an(iatom)%vxct1(:,:,:)=vxc_tmp(:,:,:)
     if (need_kxc.and.nkxc1>0) paw_an(iatom)%kxct1(:,:,:)=kxc_tmp(:,:,:)
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
&       paw_an(iatom)%lmselect,nhat1,nkxc1,mesh_size,nspden,option,&
&       pawang,pawrad(itypat),pawxcdev,rho1,0,0,vxc_tmp,xclevel,xc_denpos)
     else
       write(msg,'(2a)')ch10,' pawdenpot : Computing valence-only v_xc[n1] using angular mesh '
       call wrtout(std_out,msg,'COLL')

       call pawxc(pawtab(itypat)%coredens,eexc_val,eexcdc_val,ixc,kxc_tmp,lm_size,&
&       paw_an(iatom)%lmselect,nhat1,nkxc1,mesh_size,nspden,option,&
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
&       paw_an(iatom)%lmselect,nhat1,nkxc1,mesh_size,nspden,option,&
&       pawang,pawrad(itypat),pawxcdev,trho1,0,2*usexcnhat,vxc_tmp,xclevel,xc_denpos)
     else
       write(msg,'(2a)')ch10,' pawdenpot : Computing valence-only v_xc[tn1+nhat] using angular mesh'
       call wrtout(std_out,msg,'COLL')
       call pawxc(pawtab(itypat)%tcoredens(:,1),&
&       eexc_val,eexcdc_val,ixc,kxc_tmp,lm_size,&
&       paw_an(iatom)%lmselect,nhat1,nkxc1,mesh_size,nspden,option,&
&       pawang,pawrad(itypat),trho1,0,2*usexcnhat,vxc_tmp,xclevel,xc_denpos)
     end if
     if (option<2) then
       paw_an(iatom)%vxct1_val(:,:,:)=vxc_tmp(:,:,:)
       paw_an(iatom)%has_vxcval=2
     end if
   end if ! valence-only XC potentials

   ABI_DEALLOCATE(vxc_tmp)
   ABI_DEALLOCATE(kxc_tmp)

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
&     paw_an(iatom)%lmselect,nhat1,nkxc1,mesh_size,nspden,option,pawang,pawrad(itypat),pawxcdev,&
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
     ABI_ALLOCATE(dij_ep,(cplex*lmn2_size))
     call pawdijhartree(cplex,dij_ep,nspden,electronpositron%pawrhoij_ep(iatom),pawtab(itypat))
     if (option/=1) then
       do ispden=1,nspdiag
         jrhoij=1
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=pawrhoij(iatom)%rhoijselect(irhoij)
           ro(1)=pawrhoij(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
           electronpositron%e_paw  =electronpositron%e_paw  -ro(1)*dij_ep(klmn)
           electronpositron%e_pawdc=electronpositron%e_pawdc-ro(1)*dij_ep(klmn)
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
           klmn=pawrhoij(iatom)%rhoijselect(irhoij);kklmn=2*klmn-1
           ro(1:2)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+1,ispden)*pawtab(itypat)%dltij(klmn)
           eh2=eh2+ro(1)*paw_ij(iatom)%dijhartree(kklmn)+ro(2)*paw_ij(iatom)%dijhartree(kklmn+1)
!          Imaginary part (not used)
!          eh2=eh2+ro(2)*paw_ij(iatom)%dijhartree(kklmn)-ro(1)*paw_ij(iatom)%dijhartree(kklmn+1)
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
   keep_vhartree=(paw_an(iatom)%has_vhartree>0)
   if ((pawspnorb>0.and.ipert==0.and.ipositron/=1).or.keep_vhartree) then
     if (paw_an(iatom)%has_vhartree==0)  then
       ABI_ALLOCATE(paw_an(iatom)%vh1,(mesh_size,1,1))
     end if
     ABI_ALLOCATE(rho,(mesh_size))
     rho(1:mesh_size)=(rho1(1:mesh_size,1,1)+sqrt(four_pi)*pawtab(itypat)%coredens(1:mesh_size)) &
&     *four_pi*pawrad(itypat)%rad(1:mesh_size)**2
     call poisson(rho,0,pawrad(itypat),paw_an(iatom)%vh1(:,1,1))
     paw_an(iatom)%vh1(2:mesh_size,1,1)=(paw_an(iatom)%vh1(2:mesh_size,1,1) &
&     -sqrt(four_pi)*znucl(itypat))/pawrad(itypat)%rad(2:mesh_size)
     call pawrad_deducer0(paw_an(iatom)%vh1(:,1,1),mesh_size,pawrad(itypat))
     paw_an(iatom)%has_vhartree=2
     ABI_DEALLOCATE(rho)
   end if

!  ========= Compute PAW+U and energy contribution  =========
!  ==========================================================

   if (pawtab(itypat)%usepawu>0.and.ipert==0.and.ipositron/=1.and.option/=1.and.pawtab(itypat)%usepawu<10) then
     call pawuenergy(iatom_tot,eldaumdc,eldaumdcdc,pawprtvol,pawtab(itypat),paw_ij(iatom))
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
     if(pawrhoij(iatom)%cplex/=2) then
       msg='  pawrhoij must be complex for nuclear dipole moments !'
       MSG_BUG(msg)
     end if
     jrhoij=2 !Select imaginary part of rhoij
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
       klmn=pawrhoij(iatom)%rhoijselect(irhoij)
! select imaginary part of dijnd (only nonzero part)
       kklmn=cplex_dij*(klmn-1)+2
! because dijnd involves no spin flips, following formula is correct for all values of nspden
       enucdip=enucdip-pawrhoij(iatom)%rhoijp(jrhoij,1)*paw_ij(iatom)%dijnd(kklmn,1) &
&       *pawtab(itypat)%dltij(klmn)
       jrhoij=jrhoij+pawrhoij(iatom)%cplex
     end do
   end if

!  ========= Compute spin-orbit energy contribution  ========
!  ==========================================================

!  Compute spin-orbit contribution to Dij
   if (pawspnorb>0.and.ipert==0.and.ipositron/=1.and.&
&   (option/=2.or.pawrhoij(iatom)%cplex==2)) then
     call pawdijso(cplex_dij,paw_ij(iatom)%dijso,ndij,nspden,pawang,pawrad(itypat),pawtab(itypat), &
&     pawxcdev,spnorbscl,paw_an(iatom)%vh1,paw_an(iatom)%vxc1)
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
   if (option/=1.and.pawspnorb>0.and.ipert==0.and.ipositron/=1.and.pawrhoij(iatom)%cplex==2) then
     if(pawrhoij(iatom)%nspden/=4) then
       msg='  pawrhoij must have 4 components !'
       MSG_BUG(msg)
     end if
     jrhoij=2 !Select imaginary part of rhoij
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
       klmn=pawrhoij(iatom)%rhoijselect(irhoij)
       kklmn=cplex_dij*(klmn-1)+1
       espnorb=espnorb-pawrhoij(iatom)%rhoijp(jrhoij,3)*paw_ij(iatom)%dijso(kklmn,3) &
&       *pawtab(itypat)%dltij(klmn)
       if (cplex_dij==2) then
         kklmn=kklmn+1
         espnorb=espnorb-(pawrhoij(iatom)%rhoijp(jrhoij,2)*paw_ij(iatom)%dijso(kklmn,3) &
&         +half*pawrhoij(iatom)%rhoijp(jrhoij,4)*(paw_ij(iatom)%dijso(kklmn,1) &
&         -paw_ij(iatom)%dijso(kklmn,2))) &
&         *pawtab(itypat)%dltij(klmn)
       end if
       jrhoij=jrhoij+pawrhoij(iatom)%cplex
     end do
   end if

!  === Compute 2nd part of local exact-exchange energy and potential  ===
!  ======================================================================

   if (pawtab(itypat)%useexexch>0.and.ipert==0.and.ipositron/=1) then

     if(paw_ij(iatom)%nspden==4)  then
       msg='  Local exact-exch. not implemented for nspden=4 !'
       MSG_ERROR(msg)
     end if

     if (option<2) then
       call pawxpot(nspden,pawprtvol,pawrhoij(iatom),pawtab(itypat),paw_ij(iatom)%vpawx)
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

     if (ipositron/=1) then

!      Fock contribution to Dij
       ABI_ALLOCATE(dijfock_vv,(cplex_dij*lmn2_size,ndij))
       ABI_ALLOCATE(dijfock_cv,(cplex_dij*lmn2_size,ndij))
       call pawdijfock(cplex,cplex_dij,dijfock_vv,dijfock_cv,hyb_mixing_,hyb_mixing_sr_, &
&       ndij,nspden,nsppol,pawrhoij(iatom),pawtab(itypat))
       paw_ij(iatom)%dijfock(:,:)=dijfock_vv(:,:)+dijfock_cv(:,:)

!      Fock contribution to energy

       if (option/=1) then
         if ((cplex==1).or.(ipert==0)) then
           do ispden=1,pawrhoij(iatom)%nspden
             jrhoij=1
             do irhoij=1,pawrhoij(iatom)%nrhoijsel
               klmn=pawrhoij(iatom)%rhoijselect(irhoij)
               ro(1)=pawrhoij(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
               efockdc=efockdc+ro(1)*half*dijfock_vv(klmn,ispden)
               efock=efock+ro(1)*(half*dijfock_vv(klmn,ispden)+dijfock_cv(klmn,ispden))
               jrhoij=jrhoij+pawrhoij(iatom)%cplex
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
&               +ro(2)*(half*dijfock_vv(kklmn+1,ispden)+dijfock_cv(kklmn+1,ispden))
               jrhoij=jrhoij+pawrhoij(iatom)%cplex
             end do
           end do
         end if
       end if
!efock=zero;efockdc=zero
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
         jrhoij=jrhoij+pawrhoij(iatom)%cplex
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
     epaw=e1xc-etild1xc+eh2
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

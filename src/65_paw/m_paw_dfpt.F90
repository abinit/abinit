!!****m* m_paw_dfpt/m_paw_dfpt
!! NAME
!!  m_paw_dfpt
!!
!! FUNCTION
!!  This module contains several routines related to the 1st and 2nd order derivatives
!!    (in the DFPT approach) of PAW on-site quantities.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2020 ABINIT group (MT,AM,FJ,JWZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_dfpt

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_time, only : timab

 use defs_datatypes, only : pseudopotential_type
 use m_pawang,       only : pawang_type
 use m_pawrad,       only : pawrad_type
 use m_pawtab,       only : pawtab_type
 use m_paw_an,       only : paw_an_type
 use m_paw_ij,       only : paw_ij_type
 use m_pawcprj,      only : pawcprj_type
 use m_pawdij,       only : pawdijhartree,pawdiju_euijkl
 use m_pawrhoij,     only : pawrhoij_type,pawrhoij_free,pawrhoij_gather,pawrhoij_nullify
 use m_pawfgrtab,    only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_nullify, pawfgrtab_gather
 use m_paw_finegrid, only : pawgylm, pawrfgd_fft, pawexpiqr
 use m_pawxc,        only : pawxc_dfpt, pawxcm_dfpt
 use m_paw_denpot,   only : pawdensities,pawaccenergy,pawaccenergy_nospin
 use m_paral_atom,   only : get_my_atmtab,free_my_atmtab
 use m_atm2fft,      only : dfpt_atm2fft
 use m_distribfft,   only : distribfft_type,init_distribfft_seq,destroy_distribfft
 use m_geometry,     only : metric, stresssym
 use m_efield,       only : efield_type

 implicit none

 private

!public procedures.
 public :: pawdfptenergy ! Compute Hartree+XC PAW on-site contrib. to a 1st or 2nd-order energy
 public :: pawgrnl       ! Compute derivatives of total energy due to NL terms (PAW Dij derivatives)
 public :: dsdr_k_paw    ! Compute PAW on-site terms for forces/stresses for finite electric fields

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_dfpt/pawdfptenergy
!! NAME
!! pawdfptenergy
!!
!! FUNCTION
!! This routine compute the Hartree+XC+U PAW on-site contributions to a 1st-order or 2nd-order energy.
!!  These contributions are equal to:
!!    E_onsite=
!!       Int{ VHxc[n1_a^(j1);nc^(j1)].n1_b^(j2) }
!!      -Int{ VHxc[tild_n1_a^(j1)+hat_n1_a^(j1);tild_n_c^(j1)].(tild_n1_b+n1_b)^(j2) }
!! Some typical uses:
!!  A-Contribution to non-stationary expression of the 2nd-order total energy:
!!    In that case, n1_a^(1)[r]=n1^(j1)[r] and n1_b[r]=delta_n1^(j2)[r]
!!    where j1 and j2 are two given perturbations,
!!    and delta_n1^(j)[r] is the 1s-order density only due to change of WF overlap.
!!    See PRB 78, 035105 (2008) [[cite:Audouze2008]], Eq.(80).
!!    E_onsite=
!!       Int{ VHxc[n1^(j1);nc^(j1)].delta_n1^(j2) }
!!      -Int{ VHxc[tild_n1^(j1)+hat_n1^(j1);tild_n_c^(j1)].delta_(tild_n1+hat_n1)^(j2) }
!!  B-Contribution to first-order Fermi energy:
!!    In that case, n1_a^(1)[r]=n1^(j1)[r] and n1_b[r]=n1[r,EFermi]
!!    where j1 is the current perturbation, and n1[r,EFermi] is the density at Fermi level.
!!    E_onsite=
!!       Int{ VHxc[n1^(j1);nc^(j1)].n1[r,EFermi] }
!!      -Int{ VHxc[tild_n1^(j1)+hat_n1^(j1);tild_n_c^(j1)].(tild_n1+hat_n1)[r,EFermi] }
!!
!! INPUTS
!!  ipert1,ipert2=indexes of perturbations (j1) and (j2)
!!                if ipert2<=0, we compute a first-order energy
!!                if ipert2> 0, we compute a second-order energy
!!  ixc= choice of exchange-correlation scheme
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  ntypat=number of types of atoms in unit cell.
!!  nzlmopt_a= For the n1_a density:
!!            if -1, compute all LM-moments of the density and use non-zero LM-moments
!!            if  0, compute all LM-moments of the density and use all LM-moments
!!            if +1, compute only non-zero LM-moments of the density (stored before)
!!  nzlmopt_b= For the n1_b density:
!!            if -1, compute all LM-moments of the density and use non-zero LM-moments
!!            if  0, compute all LM-moments of the density and use all LM-moments
!!            if +1, compute only non-zero LM-moments of the density (stored before)
!!  paw_an0(natom) <type(paw_an_type)>=paw arrays for 0th-order quantities given on angular mesh
!!  paw_an1(natom) <type(paw_an_type)>=paw arrays for 1st-order quantities given on angular mesh
!!                                     This corresponds to (j1) perturbation
!!  paw_ij1(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!                                     This corresponds to (j1) perturbation
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij_a(natom) <type(pawrhoij_type)>= paw rhoij 1st-order occupancies for the (j1) perturbation
!!  pawrhoij_b(natom) <type(pawrhoij_type)>=
!!    if ipert2> 0: paw rhoij 1st-order occupancies for the (j2) perturbation corrsponding to n1_b^(j2)[r]
!!    if ipert2<=0: paw rhoij occupancies corresponding to n1_b[r]
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  xclevel= XC functional level
!!
!! OUTPUT
!!  delta_energy(2)= real and imaginary parts of contributions to non-stationary expression for the
!!              second derivative of the total energy
!!
!! SIDE EFFECTS
!!    ==== if paw_an1(:)%has_vxc<2, compute 1st-order XC potentials
!!      paw_an1(natom)%vxc1(cplex_a*mesh_size,:,nspden) =AE 1st-order XC potential Vxc^(j1)
!!      paw_an1(natom)%vxct1(cplex_a*mesh_size,:,nspden)=PS 1st-order XC potential tVxc^(j1)
!!    ==== if paw_ij1(:)%has_dijhartree<2, compute 1st-order Dij_hartree
!!      paw_ij1(natom)%dijhartree(cplex_a*lmn2_size)=Hartree contribution to Dij^(j1)
!!    ==== if paw_ij1(:)%has_dijU<2, compute 1st-order Dij_U
!!      paw_ij1(natom)%diju(cplex_a*lmn2_size)=DFT+U contribution to Dij^(j1)
!!
!! PARENTS
!!      m_dfpt_nstwf,m_dfpt_scfcv,m_dfptnl_pert
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawdfptenergy(delta_energy,ipert1,ipert2,ixc,my_natom,natom,ntypat,nzlmopt_a,nzlmopt_b,&
&                        paw_an0,paw_an1,paw_ij1,pawang,pawprtvol,pawrad,pawrhoij_a,pawrhoij_b,&
&                        pawtab,pawxcdev,xclevel, &
&                        mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ipert1,ipert2,ixc,my_natom,natom,ntypat,nzlmopt_a,nzlmopt_b
 integer,intent(in) :: pawprtvol,pawxcdev,xclevel
 integer,optional,intent(in) :: comm_atom
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(out) :: delta_energy(2)
 type(paw_an_type),intent(in) :: paw_an0(my_natom)
 type(paw_an_type),intent(inout) :: paw_an1(my_natom)
 type(paw_ij_type),intent(inout) :: paw_ij1(my_natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij_a(my_natom),pawrhoij_b(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer, parameter :: PAWU_ALGO_1=1,PAWU_ALGO_2=2
 integer :: cplex_a,cplex_b,cplex_vxc1,iatom,iatom_tot,ierr,itypat,lm_size_a,lm_size_b,mesh_size
 integer :: my_comm_atom,nspden,opt_compch,optexc,optvxc,pawu_algo,qphase_dijh1,qphase_diju1
 integer :: usecore,usepawu,usetcore,usexcnhat
 logical :: my_atmtab_allocated,non_magnetic_xc,paral_atom
 real(dp) :: compch,eexc,eexc_im
 character(len=500) :: msg
!arrays
 integer,pointer :: my_atmtab(:)
 logical,allocatable :: lmselect_a(:),lmselect_b(:),lmselect_tmp(:)
 real(dp) :: delta_energy_h(2),delta_energy_u(2),delta_energy_xc(2),tsec(2)
 real(dp),allocatable :: kxc_dum(:,:,:),nhat1(:,:,:),rho1(:,:,:),trho1(:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(567,1,tsec)

 if (.not.(ipert1==natom+1.or.ipert1==natom+10.or.ipert1==natom+11 &
& .or.ipert2==natom+1.or.ipert2==natom+10.or.ipert2==natom+11)) then
   if((abs(nzlmopt_a)/=1.and.nzlmopt_a/=0).or.(abs(nzlmopt_b)/=1.and.nzlmopt_b/=0)) then
     msg='invalid value for nzlmopt!'
     ABI_BUG(msg)
   end if
   if (my_natom>0) then
     if(paw_ij1(1)%has_dijhartree==0) then
       msg='dijhartree must be allocated!'
       ABI_BUG(msg)
     end if
     if (any(pawtab(1:ntypat)%usepawu/=0)) then
       if(paw_ij1(1)%has_dijU==0) then
         msg='dijU must be allocated!'
         ABI_BUG(msg)
       end if
     end if
     if(paw_an1(1)%has_vxc==0) then
       msg='vxc1 and vxct1 must be allocated!'
       ABI_BUG(msg)
     end if
     if(paw_an0(1)%has_kxc==0) then
       msg='kxc1 must be allocated!'
       ABI_BUG(msg)
     end if
     if ((ipert1<=natom.or.ipert1==natom+1.or.ipert1==natom+10.or.ipert1==natom+11).and.paw_an0(1)%has_kxc/=2) then
       msg='XC kernels for ground state must be in memory!'
       ABI_BUG(msg)
     end if
     if (paw_ij1(1)%qphase/=paw_an1(1)%cplex) then
       msg='paw_ij1()%qphase and paw_an1()%cplex must be equal!'
       ABI_BUG(msg)
     end if
     if (pawrhoij_a(1)%qphase<paw_an1(1)%cplex.or.pawrhoij_b(1)%qphase<paw_an1(1)%cplex) then
       msg='pawrhoij()%qphase must be >=paw_an1()%cplex!'
       ABI_BUG(msg)
     end if
     if (paw_ij1(1)%nspden/=paw_an1(1)%nspden) then
       msg='paw_ij1()%nspden and paw_an1()%nspden must be equal!'
       ABI_BUG(msg)
     end if
     if (pawrhoij_a(1)%nspden/=pawrhoij_b(1)%nspden) then
       msg='pawrhoij_a()%nspden must =pawrhoij_b()%nspden !'
       ABI_BUG(msg)
     end if
   end if
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Init contribution to 1st-order (or 2nd-order) energy
 delta_energy(1:2)=zero

!For some perturbations, nothing else to do
 if (ipert1==natom+1.or.ipert1==natom+10.or.ipert1==natom+11 .or. &
&    ipert2==natom+1.or.ipert2==natom+10.or.ipert2==natom+11) return

!Various inits
 opt_compch=0;optvxc=1;optexc=3
 usecore=0;usetcore=0  ! This is true for phonons and Efield pert.
 usexcnhat=maxval(pawtab(1:ntypat)%usexcnhat)
 delta_energy_xc(1:2)=zero
 delta_energy_h(1:2)=zero
 delta_energy_u(1:2)=zero

!================ Loop on atomic sites =======================
 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

   itypat=pawrhoij_a(iatom)%itypat
   mesh_size=pawtab(itypat)%mesh_size
   nspden=paw_an1(iatom)%nspden
   cplex_a=pawrhoij_a(iatom)%qphase
   cplex_b=pawrhoij_b(iatom)%qphase
   cplex_vxc1=paw_an1(iatom)%cplex
   qphase_dijh1=paw_ij1(iatom)%qphase
   qphase_diju1=paw_ij1(iatom)%qphase
   lm_size_a=paw_an1(iatom)%lm_size
   if (ipert2<=0) lm_size_b=paw_an0(iatom)%lm_size
   if (ipert2> 0) lm_size_b=paw_an1(iatom)%lm_size
   usepawu=pawtab(itypat)%usepawu
   pawu_algo=merge(PAWU_ALGO_1,PAWU_ALGO_2,ipert1<=0.and.ipert2<=0.and.usepawu>=0)
   non_magnetic_xc=(mod(abs(usepawu),10)==4)

!  If Vxc potentials are not in memory, compute them
   if (paw_an1(iatom)%has_vxc/=2) then
     ABI_ALLOCATE(rho1 ,(cplex_a*mesh_size,lm_size_a,nspden))
     ABI_ALLOCATE(trho1,(cplex_a*mesh_size,lm_size_a,nspden))
     ABI_ALLOCATE(nhat1,(cplex_a*mesh_size,lm_size_a,nspden*usexcnhat))
     ABI_ALLOCATE(lmselect_a,(lm_size_a))
     lmselect_a(:)=paw_an1(iatom)%lmselect(:)
     ABI_ALLOCATE(lmselect_tmp,(lm_size_a))
     lmselect_tmp(:)=.true.
     if (nzlmopt_a==1) lmselect_tmp(:)=lmselect_a(:)
!    Compute on-site 1st-order densities
     call pawdensities(compch,cplex_a,iatom_tot,lmselect_tmp,lmselect_a,&
&     lm_size_a,nhat1,nspden,nzlmopt_a,opt_compch,1-usexcnhat,-1,0,pawang,pawprtvol,&
&     pawrad(itypat),pawrhoij_a(iatom),pawtab(itypat),rho1,trho1)
     ABI_DEALLOCATE(lmselect_tmp)
!    Compute on-site 1st-order xc potentials
     if (pawxcdev/=0) then
       call pawxcm_dfpt(pawtab(itypat)%coredens,cplex_a,cplex_vxc1,eexc,ixc,paw_an0(iatom)%kxc1,&
&       lm_size_a,lmselect_a,nhat1,paw_an0(iatom)%nkxc1,non_magnetic_xc,mesh_size,nspden,optvxc,&
&       pawang,pawrad(itypat),rho1,usecore,0,&
&       paw_an1(iatom)%vxc1,xclevel)
       call pawxcm_dfpt(pawtab(itypat)%tcoredens(:,1),cplex_a,cplex_vxc1,eexc,ixc,paw_an0(iatom)%kxct1,&
&       lm_size_a,lmselect_a,nhat1,paw_an0(iatom)%nkxc1,non_magnetic_xc,mesh_size,nspden,optvxc,&
&       pawang,pawrad(itypat),trho1,usetcore,2*usexcnhat,&
&       paw_an1(iatom)%vxct1,xclevel)
     else
       call pawxc_dfpt(pawtab(itypat)%coredens,cplex_a,cplex_vxc1,eexc,ixc,paw_an0(iatom)%kxc1,&
&       lm_size_a,lmselect_a,nhat1,paw_an0(iatom)%nkxc1,non_magnetic_xc,mesh_size,nspden,optvxc,&
&       pawang,pawrad(itypat),rho1,usecore,0,&
&       paw_an0(iatom)%vxc1,paw_an1(iatom)%vxc1,xclevel)
       call pawxc_dfpt(pawtab(itypat)%tcoredens(:,1),cplex_a,cplex_vxc1,eexc,ixc,paw_an0(iatom)%kxct1,&
&       lm_size_a,lmselect_a,nhat1,paw_an0(iatom)%nkxc1,non_magnetic_xc,mesh_size,nspden,optvxc,&
&       pawang,pawrad(itypat),trho1,usetcore,2*usexcnhat,&
&       paw_an0(iatom)%vxct1,paw_an1(iatom)%vxct1,xclevel)
     end if

     paw_an1(iatom)%has_vxc=2
     ABI_DEALLOCATE(lmselect_a)
     ABI_DEALLOCATE(rho1)
     ABI_DEALLOCATE(trho1)
     ABI_DEALLOCATE(nhat1)
   end if ! has_vxc

!  Compute contribution to 1st-order (or 2nd-order) energy from 1st-order XC potential
   ABI_ALLOCATE(rho1 ,(cplex_b*mesh_size,lm_size_b,nspden))
   ABI_ALLOCATE(trho1,(cplex_b*mesh_size,lm_size_b,nspden))
   ABI_ALLOCATE(nhat1,(cplex_b*mesh_size,lm_size_b,nspden*usexcnhat))
   ABI_ALLOCATE(lmselect_b,(lm_size_b))
   if (ipert2<=0) lmselect_b(:)=paw_an0(iatom)%lmselect(:)
   if (ipert2> 0) lmselect_b(:)=paw_an1(iatom)%lmselect(:)
   ABI_ALLOCATE(lmselect_tmp,(lm_size_b))
   lmselect_tmp(:)=.true.
   if (nzlmopt_b==1) lmselect_tmp(:)=lmselect_b(:)
!  Compute on-site 1st-order densities
   call pawdensities(compch,cplex_b,iatom_tot,lmselect_tmp,lmselect_b,&
&   lm_size_b,nhat1,nspden,nzlmopt_b,opt_compch,1-usexcnhat,-1,0,pawang,pawprtvol,&
&   pawrad(itypat),pawrhoij_b(iatom),pawtab(itypat),rho1,trho1)
   ABI_DEALLOCATE(lmselect_tmp)
!  Compute contributions to 1st-order (or 2nd-order) energy
   if (pawxcdev/=0) then
     ABI_ALLOCATE(kxc_dum,(mesh_size,pawang%angl_size,0))
     call pawxcm_dfpt(pawtab(itypat)%coredens,cplex_b,cplex_vxc1,eexc,ixc,kxc_dum,&
&     lm_size_b,lmselect_b,nhat1,0,non_magnetic_xc,mesh_size,nspden,optexc,pawang,pawrad(itypat),&
&     rho1,usecore,0,paw_an1(iatom)%vxc1,xclevel,d2enxc_im=eexc_im)
     delta_energy_xc(1)=delta_energy_xc(1)+eexc
     delta_energy_xc(2)=delta_energy_xc(2)+eexc_im
     call pawxcm_dfpt(pawtab(itypat)%tcoredens(:,1),&
&     cplex_b,cplex_vxc1,eexc,ixc,kxc_dum,&
&     lm_size_b,lmselect_b,nhat1,0,non_magnetic_xc,mesh_size,nspden,optexc,pawang,pawrad(itypat),&
&     trho1,usetcore,2*usexcnhat,paw_an1(iatom)%vxct1,xclevel,&
&     d2enxc_im=eexc_im)
     ABI_DEALLOCATE(kxc_dum)
     delta_energy_xc(1)=delta_energy_xc(1)-eexc
     delta_energy_xc(2)=delta_energy_xc(2)-eexc_im
   else
     ABI_ALLOCATE(kxc_dum,(mesh_size,lm_size_b,0))
     call pawxc_dfpt(pawtab(itypat)%coredens,cplex_b,cplex_vxc1,eexc,ixc,kxc_dum,&
&     lm_size_b,lmselect_b,nhat1,0,non_magnetic_xc,mesh_size,nspden,optexc,pawang,pawrad(itypat),&
&     rho1,usecore,0,paw_an0(iatom)%vxc1,paw_an1(iatom)%vxc1,xclevel,d2enxc_im=eexc_im)
     delta_energy_xc(1)=delta_energy_xc(1)+eexc
     delta_energy_xc(2)=delta_energy_xc(2)+eexc_im
     call pawxc_dfpt(pawtab(itypat)%tcoredens(:,1),&
&     cplex_b,cplex_vxc1,eexc,ixc,kxc_dum,&
&     lm_size_b,lmselect_b,nhat1,0,non_magnetic_xc,mesh_size,nspden,optexc,pawang,pawrad(itypat),&
&     trho1,usetcore,2*usexcnhat,paw_an0(iatom)%vxct1,paw_an1(iatom)%vxct1,xclevel,&
&     d2enxc_im=eexc_im)
     ABI_DEALLOCATE(kxc_dum)
     delta_energy_xc(1)=delta_energy_xc(1)-eexc
     delta_energy_xc(2)=delta_energy_xc(2)-eexc_im
   end if
   ABI_DEALLOCATE(lmselect_b)
   ABI_DEALLOCATE(rho1)
   ABI_DEALLOCATE(trho1)
   ABI_DEALLOCATE(nhat1)

!  If Dij_hartree are not in memory, compute them
   if (paw_ij1(iatom)%has_dijhartree/=2) then
     call pawdijhartree(paw_ij1(iatom)%dijhartree,qphase_dijh1,paw_ij1(iatom)%nspden,&
&     pawrhoij_a(iatom),pawtab(itypat))
     paw_ij1(iatom)%has_dijhartree=2
   end if

!  Compute contribution to 1st-order(or 2nd-order) energy from 1st-order Hartree potential
   call pawaccenergy_nospin(delta_energy_h(1),pawrhoij_b(iatom),paw_ij1(iatom)%dijhartree, &
&                           1,qphase_dijh1,pawtab(itypat),epaw_im=delta_energy_h(2))

!  Compute contribution to 1st-order(or 2nd-order) energy from 1st-order PAW+U potential
   if (usepawu/=0.and.pawu_algo==PAWU_ALGO_2) then
!    If DijU are not in memory, compute them
     if (paw_ij1(iatom)%has_dijU/=2) then ! We force the recomputation of dijU in when cplex=2 to get diju_im
       call pawdiju_euijkl(paw_ij1(iatom)%dijU,paw_ij1(iatom)%cplex_dij,qphase_diju1,&
&                          paw_ij1(iatom)%ndij,pawrhoij_a(iatom),pawtab(itypat))
       paw_ij1(iatom)%has_dijU=2
     end if
!    Compute contribution to 1st-order(or 2nd-order) energy
     call pawaccenergy(delta_energy_u(1),pawrhoij_b(iatom),paw_ij1(iatom)%dijU,paw_ij1(iatom)%cplex_dij, &
&                      qphase_diju1,paw_ij1(iatom)%ndij,pawtab(itypat),epaw_im=delta_energy_u(2))
!    Add FLL double-counting contribution
     if (ipert1==0) then ! If j1/=0, Dij^FLL^(j1)=0 because it is constant
       call pawaccenergy_nospin(delta_energy_u(1),pawrhoij_b(iatom),pawtab(itypat)%euij_fll,1,1,&
&                               pawtab(itypat),epaw_im=delta_energy_u(2))
     end if
   end if

!  ================ End loop on atomic sites =======================
 end do

!Final building of 1st-order (or 2nd-order) energy
 delta_energy(1:2)=delta_energy_xc(1:2)+delta_energy_h(1:2)+delta_energy_u(1:2)

!Reduction in case of parallelism
 if (paral_atom) then
   call xmpi_sum(delta_energy,my_comm_atom,ierr)
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 call timab(567,2,tsec)

 DBG_EXIT("COLL")

end subroutine pawdfptenergy
!!***

!----------------------------------------------------------------------

!!****f* m_paw_dfpt/pawgrnl
!!
!! NAME
!! pawgrnl
!!
!! FUNCTION
!! PAW: Add to GRadients of total energy due to non-local term of Hamiltonian
!!      the contribution due to Dij derivatives
!! In particular, compute contribution to forces, stresses, dyn. matrix
!! Remember: Vnl=Sum_ij[|p_i>Dij<p_j|]
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dimnhat=second dimension of array nhat (0 or # of spin components)
!!  distribfft<type(distribfft_type)>=--optional-- contains all the information related
!!                                    to the FFT parallelism and plane sharing
!!  dyfr_cplex=1 if dyfrnl is real, 2 if it is complex
!!  gsqcut=Fourier cutoff on G^2 for "large sphere" of radius double that of the basis sphere
!!  mgfft=maximum size of 1D FFTs
!!  me_g0=--optional-- 1 if the current process treat the g=0 plane-wave (only needed when comm_fft is present)
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  comm_fft=--optional-- MPI communicator over FFT components (=mpi_comm_grid is not present)
!!  mpi_comm_grid=--optional-- MPI communicator over real space grid components (=comm_fft is not present)
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  nattyp(ntypat)=array describing how many atoms of each type in cell
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nhat(nfft,dimnhat)=compensation charge density on rectangular grid in real space
!!  nspden=number of spin-density components
!!  nsym=number of symmetries in space group
!!  ntypat=number of types of atoms
!!  optgr= 1 if gradients with respect to atomic position(s) have to be computed
!!  optgr2= 1 if 2nd gradients with respect to atomic position(s) have to be computed
!!  optstr= 1 if gradients with respect to strain(s) have to be computed
!!  optstr2= 1 if 2nd gradients with respect to strain(s) have to be computed
!!  paral_kgb=--optional-- 1 if "band-FFT" parallelism is activated (only needed when comm_fft is present)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(my_natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase (structure factor) information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  qphon(3)=wavevector of the phonon
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  typat(natom)=types of atoms
!!  ucvol=unit cell volume
!!  vtrial(nfft,nspden)= total local potential
!!  vxc(nfft,nspden)=XC potential
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! SIDE EFFECTS
!!  At input, this terms contain contribution from non-local projectors derivatives
!!  At output, they are updated with the contribution of Dij derivatives
!!  ==== if optgr=1 ====
!!   grnl(3*natom) =gradients of NL energy wrt atomic coordinates
!!  ==== if optstr=1 ====
!!   nlstr(6) =gradients of NL energy wrt strains
!!  ==== if optgr2=1 ====
!!   dyfrnl(dyfr_cplex,3,3,natom,natom) =2nd gradients of NL energy wrt atomic coordinates
!!  ==== if optstr=2 ====
!!    eltfrnl(6+3*natom,6)=non-symmetrized non-local contribution to the elastic tensor
!! NOTES
!!   In the case of parallelisation over atoms and calculation of dynamical matrix (optgr2=1)
!!   several data are gathered and no more distributed inside this routine.
!!
!! PARENTS
!!      m_d2frnl,m_forstr,m_scfcv_core
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawgrnl(atindx1,dimnhat,dyfrnl,dyfr_cplex,eltfrnl,grnl,gsqcut,mgfft,my_natom,natom,&
&          nattyp,nfft,ngfft,nhat,nlstr,nspden,nsym,ntypat,optgr,optgr2,optstr,optstr2,&
&          pawang,pawfgrtab,pawrhoij,pawtab,ph1d,psps,qphon,rprimd,symrec,typat,ucvol,vtrial,vxc,xred,&
&          mpi_atmtab,comm_atom,comm_fft,mpi_comm_grid,me_g0,paral_kgb,distribfft) ! optional arguments (parallelism)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dimnhat,dyfr_cplex,mgfft,my_natom,natom,nfft,nspden,nsym,ntypat
 integer,intent(in) :: optgr,optgr2,optstr,optstr2
 integer,optional,intent(in) :: me_g0,comm_atom,comm_fft,mpi_comm_grid,paral_kgb
 real(dp),intent(in) :: gsqcut,ucvol
 type(distribfft_type),optional,target,intent(in) :: distribfft
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx1(natom),nattyp(ntypat),ngfft(18)
 integer,intent(in) :: symrec(3,3,nsym),typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: nhat(nfft,dimnhat),ph1d(2,3*(2*mgfft+1)*natom),qphon(3)
 real(dp),intent(in) :: rprimd(3,3),vxc(nfft,nspden),xred(3,natom)
 real(dp),intent(in),target :: vtrial(nfft,nspden)
 real(dp),intent(inout) :: dyfrnl(dyfr_cplex,3,3,natom,natom*optgr2)
 real(dp),intent(inout) :: eltfrnl(6+3*natom,6),grnl(3*natom*optgr)
 real(dp),intent(inout) :: nlstr(6*optstr)
 type(pawfgrtab_type),target,intent(inout) :: pawfgrtab(:)
 type(pawrhoij_type),target,intent(inout) ::  pawrhoij(:)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: bufind,bufsiz,cplex,dimvtrial,eps_alpha,eps_beta,eps_gamma,eps_delta,iatm,iatom
 integer :: iatom_pawfgrtab,iatom_pawrhoij,iatom_tot,iatshft,ic,idiag,idir,ier,ilm,indx,irhoij
 integer :: isel,ishift_grhoij,ishift_gr,ishift2_gr,ishift_gr2,ishift_str,ishift_str2,ishift_str2is,ispden
 integer :: ispvtr,itypat,jatom,jatom_tot,jatm,jc,jrhoij,jtypat,klm,klmn,klmn1,ll,lm_size
 integer :: lm_sizej,lmax,lmin,lmn2_size,me_fft,mu,mua,mub,mushift,my_me_g0,my_comm_atom,my_comm_fft
 integer :: my_comm_grid,my_paral_kgb,n1,n2,n3,nfftot,nfgd,nfgd_jatom
 integer :: ngrad,ngrad_nondiag,ngradp,ngradp_nondiag,ngrhat,nsploop
 integer :: opt1,opt2,opt3,qne0,usexcnhat
 logical,parameter :: save_memory=.true.
 logical :: has_phase,my_atmtab_allocated
 logical :: paral_atom,paral_atom_pawfgrtab,paral_atom_pawrhoij,paral_grid
 real(dp) :: dlt_tmp,fact_ucvol,grhat_x,hatstr_diag,rcut_jatom,ro,ro_d,ucvol_
 character(len=500) :: msg
 type(distribfft_type),pointer :: my_distribfft
 type(pawfgrtab_type),pointer :: pawfgrtab_iatom,pawfgrtab_jatom
 type(pawrhoij_type),pointer :: pawrhoij_iatom,pawrhoij_jatom
!arrays
 integer,parameter :: alpha(9)=(/1,2,3,3,3,2,2,1,1/),beta(9)=(/1,2,3,2,1,1,3,3,2/)
 integer,parameter :: eps1(6)=(/1,2,3,2,3,1/),eps2(6)=(/1,2,3,3,1,2/)
 integer,parameter :: mu9(9)=(/1,2,3,4,5,6,4,5,6/)
 integer,allocatable :: atindx(:),atm_indx(:),mu4(:)
 integer,allocatable,target :: ifftsph_tmp(:)
 integer,ABI_CONTIGUOUS pointer :: ffti3_local(:),fftn3_distrib(:),ifft_jatom(:)
 integer, pointer :: my_atmtab(:)
 real(dp) :: gmet(3,3),gprimd(3,3),hatstr(6),rdum(1),rdum2(1),rmet(3,3),tmp(12)
 real(dp) :: work1(dyfr_cplex,3,3),work2(dyfr_cplex,3,3)
 real(dp),allocatable :: buf(:,:),buf1(:),dyfr(:,:,:,:,:),eltfr(:,:)
 real(dp),allocatable :: grhat_tmp(:,:),grhat_tmp2(:,:),hatgr(:)
 real(dp),allocatable :: prod(:,:),prodp(:,:),vloc(:),vpsp1_gr(:,:),vpsp1_str(:,:)
 real(dp),allocatable,target :: rfgd_tmp(:,:)
 real(dp),ABI_CONTIGUOUS pointer :: gylm_jatom(:,:),gylmgr_jatom(:,:,:),gylmgr2_jatom(:,:,:),expiqr_jatom(:,:)
 real(dp),ABI_CONTIGUOUS pointer :: rfgd_jatom(:,:),vtrial_(:,:)
 type(coeff2_type),allocatable :: prod_nondiag(:),prodp_nondiag(:)
 type(pawfgrtab_type),pointer :: pawfgrtab_(:),pawfgrtab_tot(:)
 type(pawrhoij_type),pointer :: pawrhoij_(:),pawrhoij_tot(:)

! *************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 qne0=0;if (qphon(1)**2+qphon(2)**2+qphon(3)**2>=1.d-15) qne0=1
 if (my_natom>0) then
   if ((optgr2==1.or.optstr2==1).and.pawrhoij(1)%ngrhoij==0) then
     msg='pawgrnl: inconsistency between variables optgr2/optstr2 and ngrhoij!'
     ABI_BUG(msg)
   end if
   if (pawfgrtab(1)%rfgd_allocated==0) then
     if ((optgr2==1.and.qne0==1).or.optstr2==1) then
       msg='pawgrnl: pawfgrtab()%rfgd array must be allocated!'
       ABI_BUG(msg)
     end if
   end if
   if (pawrhoij(1)%qphase/=1) then
     msg='pawgrnl: not supposed to be called with pawrhoij(:)%qphase=2!'
     ABI_BUG(msg)
   end if
 end if

!----------------------------------------------------------------------
!Parallelism setup

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 paral_atom_pawfgrtab=(size(pawfgrtab)/=natom)
 paral_atom_pawrhoij=(size(pawrhoij)/=natom)
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)
 if (paral_atom) then
   ABI_ALLOCATE(atm_indx,(natom))
   atm_indx=-1
   do iatom=1,my_natom
     atm_indx(my_atmtab(iatom))=iatom
   end do
 end if

!Set up parallelism over real space grid and/or FFT
 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3);nfftot=n1*n2*n3
 my_comm_grid=xmpi_comm_self;my_comm_fft=xmpi_comm_self;me_fft=0
 my_me_g0=1;my_paral_kgb=0;paral_grid=.false.;nullify(my_distribfft)
 if (present(mpi_comm_grid).or.present(comm_fft)) then
   if (present(mpi_comm_grid)) my_comm_grid=mpi_comm_grid
   if (present(comm_fft)) my_comm_fft=comm_fft
   if (.not.present(mpi_comm_grid)) my_comm_grid=comm_fft
   if (.not.present(comm_fft)) my_comm_fft=mpi_comm_grid
   paral_grid=(xmpi_comm_size(my_comm_grid)>1)
   me_fft=xmpi_comm_rank(my_comm_fft)
 end if
 if (optgr2==1.or.optstr2==1) then
   if (present(comm_fft)) then
     if ((.not.present(paral_kgb)).or.(.not.present(me_g0)).or.(.not.present(distribfft))) then
       ABI_BUG(' Need paral_kgb, me_g0 and distribfft with comm_fft !')
     end if
     my_me_g0=me_g0;my_paral_kgb=paral_kgb
     my_distribfft => distribfft
   else
     ABI_DATATYPE_ALLOCATE(my_distribfft,)
     call init_distribfft_seq(my_distribfft,'f',n2,n3,'fourdp')
   end if
   if (n2 == my_distribfft%n2_coarse) then
     fftn3_distrib => my_distribfft%tab_fftdp3_distrib
     ffti3_local => my_distribfft%tab_fftdp3_local
   else
     fftn3_distrib => my_distribfft%tab_fftdp3dg_distrib
     ffti3_local => my_distribfft%tab_fftdp3dg_local
   end if
 else
   nullify(my_distribfft,fftn3_distrib,ffti3_local)
 end if

!----------------------------------------------------------------------
!Initializations

!Compute different geometric tensors
!ucvol is not computed here but provided as input arg
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol_)
 fact_ucvol=ucvol/dble(nfftot)

!Retrieve local potential according to the use of nhat in XC
 usexcnhat=maxval(pawtab(1:ntypat)%usexcnhat)
 if (usexcnhat==0) then
   ABI_ALLOCATE(vtrial_,(nfft,1))
   dimvtrial=1
!$OMP PARALLEL DO PRIVATE(ic) SHARED(nfft,vtrial,vtrial_,vxc)
   do ic=1,nfft
     vtrial_(ic,1)=vtrial(ic,1)-vxc(ic,1)
   end do
 else
   dimvtrial=nspden
   vtrial_ => vtrial
 end if

!Initializations and allocations
 ngrhat=0;ngrad=0;ngradp=0;ngrad_nondiag=0;ngradp_nondiag=0
 ishift_grhoij=0;ishift_gr=0;ishift_gr2=0;ishift_str=0;ishift_str2=0;ishift_str2is=0;ishift2_gr=0
 cplex=1;if (qne0==1) cplex=2
 if (optgr==1) then
   ABI_ALLOCATE(hatgr,(3*natom))
   hatgr=zero
   ngrad=ngrad+3
   ngrhat=ngrhat+3
   ishift_gr2=ishift_gr2+3
 end if
 if (optgr2==1) then
   mu=min(dyfr_cplex,cplex)
   ngrad =ngrad +9
   ngradp=ngradp+3
   ngrad_nondiag = ngrad_nondiag +9*mu
   ngradp_nondiag= ngradp_nondiag+3*mu
   ngrhat= ngrhat+9*mu
 end if
 if (optstr==1) then
   hatstr=zero
   ngrad=ngrad+6
   ngrhat=ngrhat+6
   ishift_gr=ishift_gr+6
   ishift_gr2=ishift_gr2+6
   ishift_str2=ishift_str2+6
   ishift_str2is = ishift_str2is+6
 end if
 if (optstr2==1) then
   ngrad =ngrad+6*(6+3)
   ngradp=ngradp+(6+3)
   ngrad_nondiag =ngrad_nondiag+6*(6+3)
   ngradp_nondiag=ngradp_nondiag+3
   ishift2_gr=ishift2_gr+3
   ngrhat=ngrhat+6*(6+3)
   ishift_gr=ishift_gr+(6+3)
   ishift_gr2=ishift_gr2+6*(6+3)
   ishift_str2is=ishift_str2is+36
   ishift_grhoij = 6
 end if

 nsploop=nspden;if (dimvtrial<nspden) nsploop=2
 if (optgr2/=1.and.optstr2/=1) then
   ABI_ALLOCATE(grhat_tmp,(ngrhat,1))
 else
   ABI_ALLOCATE(grhat_tmp,(ngrhat,natom))
   grhat_tmp=zero
   ABI_DATATYPE_ALLOCATE(prod_nondiag,(natom))
   ABI_DATATYPE_ALLOCATE(prodp_nondiag,(natom))
   ABI_ALLOCATE(atindx,(natom))
   if(optgr2==1.or.optstr2==1)then
     ABI_ALLOCATE(vpsp1_gr,(cplex*nfft,3))
     vpsp1_gr(:,:)= zero
   end if
   if (optgr2==1) then
     ABI_ALLOCATE(dyfr,(dyfr_cplex,3,3,natom,natom))
     dyfr=zero
   end if
   if (optstr2==1) then
     ABI_ALLOCATE(vpsp1_str,(cplex*nfft,6))
     ABI_ALLOCATE(grhat_tmp2,(18,natom))
     ABI_ALLOCATE(eltfr,(6+3*natom,6))
     eltfr=zero
   end if
   ABI_ALLOCATE(mu4,(4))
   atindx(:)=0
   do iatom=1,natom
     iatm=0
     do while (atindx(iatom)==0.and.iatm<natom)
       iatm=iatm+1;if (atindx1(iatm)==iatom) atindx(iatom)=iatm
     end do
   end do
 end if

!The computation of dynamical matrix and elastic tensor requires the knowledge of
!g_l(r-R).Y_lm(r-R) and derivatives for all atoms
!Compute them here, except memory saving is activated
 if ((.not.save_memory).and.(optgr2==1.or.optstr2==1)) then
   do jatom=1,size(pawfgrtab)
     jatom_tot=jatom;if (paral_atom_pawfgrtab) jatom_tot=my_atmtab(jatom)
     pawfgrtab_jatom => pawfgrtab(jatom)
     lm_sizej=pawfgrtab_jatom%l_size**2
     opt1=0;opt2=0;opt3=0
     if (pawfgrtab_jatom%gylm_allocated==0) then
       if (allocated(pawfgrtab_jatom%gylm))  then
         ABI_DEALLOCATE(pawfgrtab_jatom%gylm)
       end if
       ABI_ALLOCATE(pawfgrtab_jatom%gylm,(pawfgrtab_jatom%nfgd,lm_sizej))
       pawfgrtab_jatom%gylm_allocated=2;opt1=1
     end if
     if (pawfgrtab_jatom%gylmgr_allocated==0) then
       if (allocated(pawfgrtab_jatom%gylmgr))  then
         ABI_DEALLOCATE(pawfgrtab_jatom%gylmgr)
       end if
       ABI_ALLOCATE(pawfgrtab_jatom%gylmgr,(3,pawfgrtab_jatom%nfgd,lm_sizej))
       pawfgrtab_jatom%gylmgr_allocated=2;opt2=1
     end if
     if (opt1+opt2+opt3>0) then
       call pawgylm(pawfgrtab_jatom%gylm,pawfgrtab_jatom%gylmgr,&
&       pawfgrtab_jatom%gylmgr2,lm_sizej,pawfgrtab_jatom%nfgd,&
&       opt1,opt2,opt3,pawtab(typat(jatom_tot)),pawfgrtab_jatom%rfgd)
     end if
     if (optgr2==1.and.qne0==1) then
       if (pawfgrtab_jatom%expiqr_allocated==0) then
         if (allocated(pawfgrtab_jatom%expiqr))  then
           ABI_DEALLOCATE(pawfgrtab_jatom%expiqr)
         end if
         pawfgrtab_jatom%expiqr_allocated=2
         ABI_ALLOCATE(pawfgrtab_jatom%expiqr,(2,nfgd))
         call pawexpiqr(pawfgrtab_jatom%expiqr,gprimd,pawfgrtab_jatom%nfgd,&
&         qphon,pawfgrtab_jatom%rfgd,xred(:,jatom_tot))
       end if
     end if
   end do
 end if

!The computation of dynamical matrix and elastic tensor might require some communications
 if ((optgr2==1.or.optstr2==1).and.paral_atom.and.paral_atom_pawfgrtab.and.(.not.save_memory)) then
   ABI_DATATYPE_ALLOCATE(pawfgrtab_tot,(natom))
   call pawfgrtab_nullify(pawfgrtab_tot)
   call pawfgrtab_gather(pawfgrtab,pawfgrtab_tot,my_comm_atom,ier,mpi_atmtab=my_atmtab)
 else
   pawfgrtab_tot => pawfgrtab
 end if
 if ((optgr2==1.or.optstr2==1).and.paral_atom.and.paral_atom_pawrhoij) then
   ABI_DATATYPE_ALLOCATE(pawrhoij_tot,(natom))
   call pawrhoij_nullify(pawrhoij_tot)
   call pawrhoij_gather(pawrhoij,pawrhoij_tot,-1,my_comm_atom, &
&   with_rhoijres=.false.,with_rhoij_=.false.,with_lmnmix=.false.)
 else
   pawrhoij_tot => pawrhoij
 end if

 if (save_memory) then
   pawfgrtab_ => pawfgrtab
   pawrhoij_  => pawrhoij
 else
   pawfgrtab_ => pawfgrtab_tot
   pawrhoij_  => pawrhoij_tot
 end if

!----------------------------------------------------------------------
!Loops over types and atoms

 iatshft=0
 do itypat=1,ntypat

   lmn2_size=pawtab(itypat)%lmn2_size
   lm_size=pawtab(itypat)%lcut_size**2

   do iatm=iatshft+1,iatshft+nattyp(itypat)

     iatom_tot=atindx1(iatm)
     iatom=iatom_tot
     if (paral_atom) then
       if (save_memory.or.(optgr2/=1.and.optstr2/=1)) iatom=atm_indx(iatom_tot)
     end if

     if (iatom==-1) cycle
     iatom_pawfgrtab=iatom_tot;if (paral_atom_pawfgrtab) iatom_pawfgrtab=iatom
     iatom_pawrhoij =iatom_tot;if (paral_atom_pawrhoij)  iatom_pawrhoij =iatom
     pawfgrtab_iatom => pawfgrtab_(iatom_pawfgrtab)
     pawrhoij_iatom  => pawrhoij_(iatom_pawrhoij)

     idiag=1;if (optgr2==1.or.optstr2==1) idiag=iatm
     nfgd=pawfgrtab_iatom%nfgd

     ABI_ALLOCATE(vloc,(nfgd))
     if (ngrad>0)  then
       ABI_ALLOCATE(prod,(ngrad,lm_size))
     end if
     if (ngradp>0)  then
       ABI_ALLOCATE(prodp,(ngradp,lm_size))
     end if
     if (ngrad_nondiag>0.and.ngradp_nondiag>0) then
       do jatm=1,natom
         jtypat=typat(atindx1(jatm))
         lm_sizej=pawtab(jtypat)%lcut_size**2
         ABI_ALLOCATE(prod_nondiag(jatm)%value,(ngrad_nondiag,lm_sizej))
         ABI_ALLOCATE(prodp_nondiag(jatm)%value,(ngradp_nondiag,lm_sizej))
       end do
     end if

     grhat_tmp=zero
     if(optstr2==1) grhat_tmp2=zero

!    ------------------------------------------------------------------
!    Compute some useful data

!    Eventually compute g_l(r).Y_lm(r) derivatives for the current atom (if not already done)
     if ((optgr==1.or.optstr==1).and.(optgr2/=1).and.(optstr2/=1)) then
       if (pawfgrtab_iatom%gylmgr_allocated==0) then
         if (allocated(pawfgrtab_iatom%gylmgr))  then
           ABI_DEALLOCATE(pawfgrtab_iatom%gylmgr)
         end if
         ABI_ALLOCATE(pawfgrtab_iatom%gylmgr,(3,pawfgrtab_iatom%nfgd,lm_size))
         pawfgrtab_iatom%gylmgr_allocated=2
         call pawgylm(rdum,pawfgrtab_iatom%gylmgr,rdum2,lm_size,pawfgrtab_iatom%nfgd,&
&         0,1,0,pawtab(itypat),pawfgrtab_iatom%rfgd)
       end if

     end if
     if (optgr2==1.or.optstr2==1) then
       opt1=0;opt2=0;opt3=0
       if (pawfgrtab_iatom%gylm_allocated==0) then
         if (allocated(pawfgrtab_iatom%gylm))  then
           ABI_DEALLOCATE(pawfgrtab_iatom%gylm)
         end if
         ABI_ALLOCATE(pawfgrtab_iatom%gylm,(pawfgrtab_iatom%nfgd,lm_size))
         pawfgrtab_iatom%gylm_allocated=2;opt1=1
       end if
       if (pawfgrtab_iatom%gylmgr_allocated==0) then
         if (allocated(pawfgrtab_iatom%gylmgr))  then
           ABI_DEALLOCATE(pawfgrtab_iatom%gylmgr)
         end if
         ABI_ALLOCATE(pawfgrtab_iatom%gylmgr,(3,pawfgrtab_iatom%nfgd,lm_size))
         pawfgrtab_iatom%gylmgr_allocated=2;opt2=1
       end if
       if (pawfgrtab_iatom%gylmgr2_allocated==0) then
         if (allocated(pawfgrtab_iatom%gylmgr2))  then
           ABI_DEALLOCATE(pawfgrtab_iatom%gylmgr2)
         end if
         ABI_ALLOCATE(pawfgrtab_iatom%gylmgr2,(6,pawfgrtab_iatom%nfgd,lm_size))
         pawfgrtab_iatom%gylmgr2_allocated=2;opt3=1
       end if
       if (opt1+opt2+opt3>0) then
         call pawgylm(pawfgrtab_iatom%gylm,pawfgrtab_iatom%gylmgr,&
&         pawfgrtab_iatom%gylmgr2,lm_size,pawfgrtab_iatom%nfgd,&
&         opt1,opt2,opt3,pawtab(itypat),pawfgrtab_iatom%rfgd)
       end if
     end if

!    Eventually compute exp(-i.q.r) factors for the current atom (if not already done)
     if (optgr2==1.and.qne0==1.and.(pawfgrtab_iatom%expiqr_allocated==0)) then
       if (allocated(pawfgrtab_iatom%expiqr))  then
         ABI_DEALLOCATE(pawfgrtab_iatom%expiqr)
       end if
       ABI_ALLOCATE(pawfgrtab_iatom%expiqr,(2,nfgd))
       call pawexpiqr(pawfgrtab_iatom%expiqr,gprimd,nfgd,qphon,&
&       pawfgrtab_iatom%rfgd,xred(:,iatom))
       pawfgrtab_iatom%expiqr_allocated=2
     end if
     has_phase=(optgr2==1.and.pawfgrtab_iatom%expiqr_allocated/=0)

!    Eventually compute 1st-order potential
     if (optgr2==1.or.optstr2==1) then
       call dfpt_atm2fft(atindx,cplex,gmet,gprimd,gsqcut,idir,iatom_tot,&
&       mgfft,psps%mqgrid_vl,natom,3,nfft,ngfft,ntypat,ph1d,&
&       psps%qgrid_vl,qphon,typat,ucvol,psps%usepaw,xred,psps,pawtab,atmvlocr1=vpsp1_gr,&
&       vspl=psps%vlspl,comm_fft=my_comm_fft,me_g0=my_me_g0,&
&       paral_kgb=my_paral_kgb,distribfft=my_distribfft)
       if (cplex==1) then
         do ic=1,nfft
           tmp(1:3)=vpsp1_gr(ic,1:3)
           do mu=1,3
             vpsp1_gr(ic,mu)=-(gprimd(mu,1)*tmp(1)+gprimd(mu,2)*tmp(2)+gprimd(mu,3)*tmp(3))
           end do
         end do
       else ! cplex=2
         do ic=1,nfft
           jc=2*ic;tmp(1:3)=vpsp1_gr(jc-1,1:3);tmp(4:6)=vpsp1_gr(jc,1:3)
           do mu=1,3
             vpsp1_gr(jc-1,mu)=-(gprimd(mu,1)*tmp(1)+gprimd(mu,2)*tmp(2)+gprimd(mu,3)*tmp(3))
             vpsp1_gr(jc  ,mu)=-(gprimd(mu,1)*tmp(4)+gprimd(mu,2)*tmp(5)+gprimd(mu,3)*tmp(6))
           end do
         end do
       end if
     end if
     if (optstr2==1) then
       vpsp1_str(:,:) = zero
       call dfpt_atm2fft(atindx,cplex,gmet,gprimd,gsqcut,idir,natom+3,&
&       mgfft,psps%mqgrid_vl,natom,6,nfft,ngfft,ntypat,&
&       ph1d,psps%qgrid_vl,qphon,typat,ucvol,psps%usepaw,xred,psps,pawtab,atmvlocr1=vpsp1_str,&
&       vspl=psps%vlspl,comm_fft=my_comm_fft,me_g0=my_me_g0,&
&       paral_kgb=my_paral_kgb,distribfft=my_distribfft)
     end if

!    ------------------------------------------------------------------
!    Loop over spin components

     do ispden=1,nsploop

!      ----- Retrieve potential (subtle if nspden=4 ;-)
       if (nspden/=4) then
         ispvtr=min(dimvtrial,ispden)
         do ic=1,nfgd
           jc = pawfgrtab_iatom%ifftsph(ic)
           vloc(ic)=vtrial_(jc,ispvtr)
         end do
       else
         if (ispden==1) then
           ispvtr=min(dimvtrial,2)
           do ic=1,nfgd
             jc=pawfgrtab_iatom%ifftsph(ic)
             vloc(ic)=half*(vtrial_(jc,1)+vtrial_(jc,ispvtr))
           end do
         else if (ispden==4) then
           ispvtr=min(dimvtrial,2)
           do ic=1,nfgd
             jc=pawfgrtab_iatom%ifftsph(ic)
             vloc(ic)=half*(vtrial_(jc,1)-vtrial_(jc,ispvtr))
           end do
         else if (ispden==2) then
           ispvtr=min(dimvtrial,3)
           do ic=1,nfgd
             jc=pawfgrtab_iatom%ifftsph(ic)
             vloc(ic)=vtrial_(jc,ispvtr)
           end do
         else ! ispden=3
           ispvtr=min(dimvtrial,4)
           do ic=1,nfgd
             jc=pawfgrtab_iatom%ifftsph(ic)
             vloc(ic)=-vtrial_(jc,ispvtr)
           end do
         end if
       end if

!      -----------------------------------------------------------------------
!      ----- Compute projected scalars (integrals of vloc and Q_ij^hat) ------
!      ----- and/or their derivatives ----------------------------------------

       if (ngrad>0) prod=zero
       if (ngradp>0) prodp=zero

!      ==== Contribution to forces ====
       if (optgr==1) then
         do ilm=1,lm_size
           do ic=1,pawfgrtab_iatom%nfgd
             do mu=1,3
               prod(mu+ishift_gr,ilm)=prod(mu+ishift_gr,ilm)-&
&               vloc(ic)*pawfgrtab_iatom%gylmgr(mu,ic,ilm)
             end do
           end do
         end do
       end if ! optgr

!      ==== Contribution to stresses ====
       if (optstr==1) then
         do ilm=1,lm_size
           do ic=1,pawfgrtab_iatom%nfgd
             jc=pawfgrtab_iatom%ifftsph(ic)
             do mu=1,6
               mua=alpha(mu);mub=beta(mu)
               prod(mu+ishift_str,ilm)=prod(mu+ishift_str,ilm) &
&               +half*vloc(ic)&
&               *(pawfgrtab_iatom%gylmgr(mua,ic,ilm)*pawfgrtab_iatom%rfgd(mub,ic)&
&               +pawfgrtab_iatom%gylmgr(mub,ic,ilm)*pawfgrtab_iatom%rfgd(mua,ic))
             end do
           end do
         end do
       end if ! optstr

!      ==== Diagonal contribution to frozen wf part of dyn. matrix ====
       if (optgr2==1) then
!        Diagonal contribution
         do ilm=1,lm_size
           do ic=1,pawfgrtab_iatom%nfgd
             do mu=1,9
               prod(ishift_gr2+mu,ilm)=prod(ishift_gr2+mu,ilm) &
&               +half*vloc(ic)*pawfgrtab_iatom%gylmgr2(mu9(mu),ic,ilm)
             end do
             do mu=1,3
               prodp(ishift_gr+mu,ilm)=prodp(ishift_gr+mu,ilm) &
&               -vloc(ic)*pawfgrtab_iatom%gylmgr(mu,ic,ilm)
             end do
           end do
         end do
       end if ! optgr2

!      ==== Diagonal contribution to elastic tensor ====
       if (optstr2==1) then
         do ilm=1,lm_size
           do ic=1,pawfgrtab_iatom%nfgd
             mu=1
             jc=pawfgrtab_iatom%ifftsph(ic)
             do mua=1,6
               eps_alpha=eps1(mua);eps_beta=eps2(mua);
               do mub=1,6
                 eps_gamma=eps1(mub);eps_delta=eps2(mub);
                 mu4 = zero
                 call pawgrnl_convert(mu4,eps_alpha,eps_beta,eps_gamma,eps_delta)
!                v_loc*d2glylm
                 prod(ishift_str2+mu,ilm)=prod(ishift_str2+mu,ilm) + half*half*vloc(ic)*( &
&                 pawfgrtab_iatom%rfgd(eps_beta,ic)*pawfgrtab_iatom%rfgd(eps_gamma,ic)*&
&                 pawfgrtab_iatom%gylmgr2(mu4(2),ic,ilm)&
&                 +pawfgrtab_iatom%rfgd(eps_alpha,ic)*pawfgrtab_iatom%rfgd(eps_gamma,ic)*&
                 pawfgrtab_iatom%gylmgr2(mu4(4),ic,ilm)&
&                 +pawfgrtab_iatom%rfgd(eps_beta,ic) *pawfgrtab_iatom%rfgd(eps_delta,ic)*&
                 pawfgrtab_iatom%gylmgr2(mu4(1),ic,ilm)&
&                 +pawfgrtab_iatom%rfgd(eps_alpha,ic)*pawfgrtab_iatom%rfgd(eps_delta,ic)*&
                 pawfgrtab_iatom%gylmgr2(mu4(3),ic,ilm))
                 if(eps_gamma==eps_beta)then
                   prod(ishift_str2+mu,ilm)=prod(ishift_str2+mu,ilm) &
&                   +half*half*vloc(ic)*(pawfgrtab_iatom%gylmgr(eps_delta,ic,ilm)*pawfgrtab_iatom%rfgd(eps_alpha,ic))
                 end if
                 if(eps_gamma==eps_alpha)then
                   prod(ishift_str2+mu,ilm)=prod(ishift_str2+mu,ilm) &
&                   +half*half*vloc(ic)*(pawfgrtab_iatom%gylmgr(eps_delta,ic,ilm)*pawfgrtab_iatom%rfgd(eps_beta,ic))
                 end if
                 if(eps_delta==eps_beta)then
                   prod(ishift_str2+mu,ilm)=prod(ishift_str2+mu,ilm) &
&                   +half*half*vloc(ic)*(pawfgrtab_iatom%gylmgr(eps_gamma,ic,ilm)*pawfgrtab_iatom%rfgd(eps_alpha,ic))
                 end if
                 if(eps_delta==eps_alpha)then
                   prod(ishift_str2+mu,ilm)=prod(ishift_str2+mu,ilm) &
&                   +half*half*vloc(ic)*(pawfgrtab_iatom%gylmgr(eps_gamma,ic,ilm)*pawfgrtab_iatom%rfgd(eps_beta,ic))
                 end if
!                d(vloc)/d(eps_gammadelta) * d(gylm)/d(eps_alphabeta)
                 prod(ishift_str2+mu,ilm)=prod(ishift_str2+mu,ilm)&
&                 +vpsp1_str(jc,mub)*half*(&
&                 (pawfgrtab_iatom%gylmgr(eps_alpha,ic,ilm)*pawfgrtab_iatom%rfgd(eps_beta,ic)&
&                 +pawfgrtab_iatom%gylmgr(eps_beta,ic,ilm) *pawfgrtab_iatom%rfgd(eps_alpha,ic)))
!                d(vloc)/d(eps_alphabeta)  * d(gylm)/d(eps_gammadelta)
                 prod(ishift_str2+mu,ilm)=prod(ishift_str2+mu,ilm)&
&                 +vpsp1_str(jc,mua)*half*(&
&                 (pawfgrtab_iatom%gylmgr(eps_gamma,ic,ilm)*pawfgrtab_iatom%rfgd(eps_delta,ic)&
&                 +pawfgrtab_iatom%gylmgr(eps_delta,ic,ilm)*pawfgrtab_iatom%rfgd(eps_gamma,ic)))
!                delta_alphabeta * dv_loc/depsgammadelta * (gylm)
                 if (mua<=3) then
                   prod(ishift_str2+mu,ilm)=prod(ishift_str2+mu,ilm) &
&                   +vpsp1_str(jc,mub)*pawfgrtab_iatom%gylm(ic,ilm)
                 end if
!                delta_gammadelta * dv_loc/depsalphabeta * (gylm)
                 if (mub<=3) then
                   prod(ishift_str2+mu,ilm)=prod(ishift_str2+mu,ilm) &
&                   +vpsp1_str(jc,mua) * pawfgrtab_iatom%gylm(ic,ilm)
                 end if
!                delta_gammadelta * v_loc * d(gylm)/d(eps_alphabeta)
                 if (mub<=3) then
                   prod(ishift_str2+mu,ilm)=prod(ishift_str2+mu,ilm) &
&                   +half*vloc(ic)&
&                   *(pawfgrtab_iatom%gylmgr(eps_beta,ic,ilm)*pawfgrtab_iatom%rfgd(eps_alpha,ic)&
&                   + pawfgrtab_iatom%gylmgr(eps_alpha,ic,ilm)*pawfgrtab_iatom%rfgd(eps_beta,ic))
                 end if
!                delta_alphabeta * v_loc * d(gylm)/d(eps_gammadelta)
                 if (mua<=3) then
                   prod(ishift_str2+mu,ilm)=prod(ishift_str2+mu,ilm) &
&                   +half*vloc(ic)&
&                   *(pawfgrtab_iatom%gylmgr(eps_gamma,ic,ilm)*pawfgrtab_iatom%rfgd(eps_delta,ic)&
&                   + pawfgrtab_iatom%gylmgr(eps_delta,ic,ilm)*pawfgrtab_iatom%rfgd(eps_gamma,ic))
                 end if
!                delta_gammadelta delta_alphabeta * v_loc * (gylm)
                 if (mua<=3.and.mub<=3) then
                   prod(ishift_str2+mu,ilm)=prod(ishift_str2+mu,ilm) &
&                   +vloc(ic)*pawfgrtab_iatom%gylm(ic,ilm)
                 end if
                 mu=mu+1
               end do !end loop mub
             end do !end loop mua
!            vloc * d(gylm)/d(eps_alphabeta)
             do mu=1,6
               mua=alpha(mu);mub=beta(mu)
               prodp(ishift_str2+mu,ilm)=prodp(ishift_str2+mu,ilm)&
&               +half*vloc(ic)*&
&               (pawfgrtab_iatom%gylmgr(mua,ic,ilm)*pawfgrtab_iatom%rfgd(mub,ic)&
&               +pawfgrtab_iatom%gylmgr(mub,ic,ilm)*pawfgrtab_iatom%rfgd(mua,ic))
!              d(vloc)/d(eps_alphabeta or gammadelta) * gylm
               prodp(ishift_str2+mu,ilm)=prodp(ishift_str2+mu,ilm)&
&               +vpsp1_str(jc,mu)*pawfgrtab_iatom%gylm(ic,ilm)
!              delta_alphabeta * vloc * gylm
               if (mu<=3) then
                 prodp(ishift_str2+mu,ilm)=prodp(ishift_str2+mu,ilm)&
&                 +vloc(ic)*pawfgrtab_iatom%gylm(ic,ilm)
               end if

!              INTERNAL STRAIN CONTRIBUTION:
               do idir=1,3
!                v_loc*d2glylm/dR contribution:
                 eps_alpha=alpha(mu);eps_beta=beta(mu);
                 call pawgrnl_convert(mu4,eps_alpha,eps_beta,idir,idir)
                 prod(ishift_str2is+(mu-1)*3+idir,ilm)=prod(ishift_str2is+(mu-1)*3+idir,ilm)&
&                 -half*vloc(ic)&
&                 *(pawfgrtab_iatom%gylmgr2(mu4(3),ic,ilm)*pawfgrtab_iatom%rfgd(eps_alpha,ic)&
&                 +pawfgrtab_iatom%gylmgr2(mu4(1),ic,ilm)*pawfgrtab_iatom%rfgd(eps_beta,ic))
                 if (idir==eps_beta)then
                   prod(ishift_str2is+(mu-1)*3+idir,ilm)=prod(ishift_str2is+(mu-1)*3+idir,ilm)&
&                   -half*vloc(ic)*(pawfgrtab_iatom%gylmgr(eps_alpha,ic,ilm))
                 end if
                 if (idir==eps_alpha)then
                   prod(ishift_str2is+(mu-1)*3+idir,ilm)=prod(ishift_str2is+(mu-1)*3+idir,ilm)&
&                   -half*vloc(ic)*(pawfgrtab_iatom%gylmgr(eps_beta,ic,ilm))
                 end if
!                delta_gammadelta * v_loc * d(gylm)/dR
                 if (mu<=3) then
                   prod(ishift_str2is+(mu-1)*3+idir,ilm)=prod(ishift_str2is+(mu-1)*3+idir,ilm)-&
                   vloc(ic)*pawfgrtab_iatom%gylmgr(idir,ic,ilm)
                 end if
!                dv_loc/deps_alph_beta * d(gylm)/dR
                 prod(ishift_str2is+(mu-1)*3+idir,ilm)=prod(ishift_str2is+(mu-1)*3+idir,ilm)-&
                 vpsp1_str(jc,mu)*pawfgrtab_iatom%gylmgr(idir,ic,ilm)
               end do
             end do
             do idir=1,3
!              v_loc * d(gylm)/dR
               prodp(6+idir,ilm) = prodp(6+idir,ilm)-vloc(ic)*pawfgrtab_iatom%gylmgr(idir,ic,ilm)
             end do !end loop idir
!            END INTERNAL STRAIN CONTRIBUTION

           end do
         end do
       end if !optstr2

!      Off-diagonal contributions
       if (optgr2==1.or.optstr2==1) then
         do jatm=1,natom
           jatom_tot=atindx1(jatm);jtypat=typat(jatom_tot)
           jatom=jatom_tot;if (paral_atom.and.save_memory) jatom=atm_indx(jatom_tot)
           lm_sizej=pawtab(jtypat)%lcut_size**2

!          Retrieve data for the atom j
           if (save_memory.and.jatom/=iatom) then
             rcut_jatom=pawtab(jtypat)%rshp
             call pawrfgd_fft(ifftsph_tmp,gmet,n1,n2,n3,nfgd_jatom,rcut_jatom,rfgd_tmp,rprimd,&
&             ucvol,xred(:,jatom_tot),fft_distrib=fftn3_distrib,fft_index=ffti3_local,me_fft=me_fft)
             ifft_jatom => ifftsph_tmp ; rfgd_jatom => rfgd_tmp
             ABI_ALLOCATE(gylm_jatom,(nfgd_jatom,lm_sizej))
             ABI_ALLOCATE(gylmgr_jatom,(3,nfgd_jatom,lm_sizej))
             opt1=1;opt2=1;opt3=0;gylmgr2_jatom=>gylmgr_jatom
             call pawgylm(gylm_jatom,gylmgr_jatom,gylmgr2_jatom,lm_sizej,nfgd_jatom,&
&             opt1,opt2,opt3,pawtab(typat(jatom_tot)),rfgd_jatom)
             if (optgr2==1.and.qne0==1) then
               ABI_ALLOCATE(expiqr_jatom,(2,nfgd_jatom))
               call pawexpiqr(expiqr_jatom,gprimd,nfgd_jatom,qphon,rfgd_jatom,xred(:,jatom_tot))
             end if
           else
             pawfgrtab_jatom => pawfgrtab_tot(jatom)
             nfgd_jatom      =  pawfgrtab_jatom%nfgd
             ifft_jatom      => pawfgrtab_jatom%ifftsph
             rfgd_jatom      => pawfgrtab_jatom%rfgd
             gylm_jatom      => pawfgrtab_jatom%gylm
             gylmgr_jatom    => pawfgrtab_jatom%gylmgr
             gylmgr2_jatom   => pawfgrtab_jatom%gylmgr2
             expiqr_jatom    => pawfgrtab_jatom%expiqr
           end if

!          ==== Off-diagonal contribution to frozen wf part of dyn. matrix ====
           if (optgr2==1) then
             mu = min(dyfr_cplex,cplex)
             prod_nondiag(jatm)%value(ishift_gr2+1:ishift_gr2+(9*mu),:) = zero
             prodp_nondiag(jatm)%value(ishift2_gr+1:ishift2_gr +(3*mu),:) = zero
             if (has_phase.or.cplex==2) then
               if (dyfr_cplex==1.or.cplex==1) then
                 do ilm=1,lm_sizej
                   do ic=1,nfgd_jatom
                     jc=2*ifft_jatom(ic)
                     tmp(1:3)=vpsp1_gr(jc-1,1:3)*expiqr_jatom(1,ic) &
&                     -vpsp1_gr(jc  ,1:3)*expiqr_jatom(2,ic)
                     do mu=1,9
                       mua=alpha(mu);mub=beta(mu)
                       prod_nondiag(jatm)%value(ishift_gr2+mu,ilm)=prod_nondiag(jatm)%value(ishift_gr2+mu,ilm) &
&                       +tmp(mua)*gylmgr_jatom(mub,ic,ilm)
                     end do
                     do mu=1,3
                       prodp_nondiag(jatm)%value(ishift2_gr+mu,ilm)=prodp_nondiag(jatm)%value(ishift2_gr+mu,ilm) &
&                       -tmp(mu)*gylm_jatom(ic,ilm)
                     end do
                   end do
                 end do
               else
                 do ilm=1,lm_sizej
                   do ic=1,nfgd_jatom
                     jc=2*ifft_jatom(ic)
                     tmp(1:3)=vpsp1_gr(jc-1,1:3)*expiqr_jatom(1,ic) &
&                     -vpsp1_gr(jc  ,1:3)*expiqr_jatom(2,ic)
                     tmp(4:6)=vpsp1_gr(jc-1,1:3)*expiqr_jatom(2,ic) &
&                     +vpsp1_gr(jc  ,1:3)*expiqr_jatom(1,ic)
                     do mu=1,9
                       mua=alpha(mu);mub=beta(mu)
                       prod_nondiag(jatm)%value(ishift_gr2+mu,ilm)=prod_nondiag(jatm)%value(ishift_gr2+mu,ilm) &
&                       +tmp(mua  )*gylmgr_jatom(mub,ic,ilm)
                       prod_nondiag(jatm)%value(ishift_gr2+9+mu,ilm)=prod_nondiag(jatm)%value(ishift_gr2+9+mu,ilm) &
&                       +tmp(3+mua)*gylmgr_jatom(mub,ic,ilm)
                     end do
                     do mu=1,3
                       prodp_nondiag(jatm)%value(ishift2_gr+mu,ilm)=prodp_nondiag(jatm)%value(ishift2_gr+mu,ilm) &
&                       -tmp(  mu)*gylm_jatom(ic,ilm)
                       prodp_nondiag(jatm)%value(ishift2_gr+3+mu,ilm)=prodp_nondiag(jatm)%value(ishift2_gr+3+mu,ilm) &
&                       -tmp(3+mu)*gylm_jatom(ic,ilm)
                     end do
                   end do
                 end do
               end if
             else ! no phase
               do ilm=1,lm_sizej
                 do ic=1,nfgd_jatom
                   jc=ifft_jatom(ic)
                   do mu=1,9
                     mua=alpha(mu);mub=beta(mu)
                     prod_nondiag(jatm)%value(ishift_gr2+mu,ilm)=prod_nondiag(jatm)%value(ishift_gr2+mu,ilm) &
&                     +vpsp1_gr(jc,mua)*gylmgr_jatom(mub,ic,ilm)
                   end do
                   do mu=1,3
                     prodp_nondiag(jatm)%value(ishift2_gr+mu,ilm)=prodp_nondiag(jatm)%value(ishift2_gr+mu,ilm) &
&                     -vpsp1_gr(jc,mu)*gylm_jatom(ic,ilm)
                   end do
                 end do
               end do
             end if
           end if ! optgr2

!          ==== Off-diagonal contribution to elastic tensor ====
           if (optstr2==1) then
             prod_nondiag(jatm)%value(ishift_str2is+1:ishift_str2is+18,:)=zero
             prodp_nondiag(jatm)%value(1:3,:)=zero
             do ilm=1,lm_sizej
               do ic=1,nfgd_jatom
                 mu=1;jc=ifft_jatom(ic)
!                INTERNAL STRAIN CONTRIBUTION:
                 do mua=1,6
                   eps_alpha=eps1(mua);eps_beta=eps2(mua);
!                  d(-vloc)/dR * d(gylm)/d(eps_gamma_delta)
                   do idir=1,3
                     prod_nondiag(jatm)%value(ishift_str2is+(mua-1)*3+idir,ilm)=&
&                     prod_nondiag(jatm)%value(ishift_str2is+(mua-1)*3+idir,ilm)&
&                     -vpsp1_gr(jc,idir)*half&
&                     *(gylmgr_jatom(eps_alpha,ic,ilm) * rfgd_jatom(eps_beta ,ic)&
&                     +gylmgr_jatom(eps_beta ,ic,ilm) * rfgd_jatom(eps_alpha,ic))
!                    delta_alphabeta * d(-v_loc/dr) * gylm
                     if (mua<=3) then
                       prod_nondiag(jatm)%value(ishift_str2is+(mua-1)*3+idir,ilm)=&
&                       prod_nondiag(jatm)%value(ishift_str2is+(mua-1)*3+idir,ilm)&
&                       -vpsp1_gr(jc,idir)*gylm_jatom(ic,ilm)
                     end if
                   end do ! dir
                 end do ! mua
                 do idir=1,3
!                  d(-v_loc/dr) * gylm
                   prodp_nondiag(jatm)%value(idir,ilm) = prodp_nondiag(jatm)%value(idir,ilm)&
&                   -vpsp1_gr(jc,idir)*gylm_jatom(ic,ilm)
                 end do !end loop idir
!                END INTERNAL STRAIN CONTRIBUTION
               end do
             end do
           end if ! optstr2

!          Release temp memory allocated for atom j
           if (save_memory.and.jatom/=iatom) then
             ABI_DEALLOCATE(ifftsph_tmp)
             ABI_DEALLOCATE(rfgd_tmp)
             ABI_DEALLOCATE(gylm_jatom)
             ABI_DEALLOCATE(gylmgr_jatom)
             if (optgr2==1.and.qne0==1) then
               ABI_DEALLOCATE(expiqr_jatom)
             end if
           end if

         end do ! loop on atoms j
       end if ! optgr2 or optstr2

!      --- Apply scaling factor on integrals ---
       if (ngrad >0) prod (:,:)=prod (:,:)*fact_ucvol
       if (ngradp>0) prodp(:,:)=prodp(:,:)*fact_ucvol
       if (ngrad_nondiag>0) then
         do jatm=1,natom
           prod_nondiag(jatm)%value(:,:)=prod_nondiag(jatm)%value(:,:)*fact_ucvol
         end do
       end if
       if (ngradp_nondiag>0) then
         do jatm=1,natom
           prodp_nondiag(jatm)%value(:,:)=prodp_nondiag(jatm)%value(:,:)*fact_ucvol
         end do
       end if

!      --- Reduction in case of parallelization ---
       if (paral_grid) then
         if (ngrad>0) then
           call xmpi_sum(prod,my_comm_grid,ier)
         end if
         if (ngradp>0) then
           call xmpi_sum(prodp,my_comm_grid,ier)
         end if
         if (ngrad_nondiag>0.or.ngradp_nondiag>0) then
           bufsiz=0;bufind=0
           do jatm=1,natom
             jtypat=typat(atindx1(jatm))
             bufsiz=bufsiz+pawtab(jtypat)%lcut_size**2
           end do
           ABI_ALLOCATE(buf,(ngrad_nondiag+ngradp_nondiag,bufsiz))
           do jatm=1,natom
             jtypat=typat(atindx1(jatm))
             lm_sizej=pawtab(jtypat)%lcut_size**2
             if (ngrad_nondiag> 0) buf(1:ngrad_nondiag,bufind+1:bufind+lm_sizej)= &
&             prod_nondiag(jatm)%value(:,:)
             if (ngradp_nondiag>0) buf(ngrad_nondiag+1:ngrad_nondiag+ngradp_nondiag, &
&             bufind+1:bufind+lm_sizej)=prodp_nondiag(jatm)%value(:,:)
             bufind=bufind+lm_sizej*(ngrad_nondiag+ngradp_nondiag)
           end do
           call xmpi_sum(buf,my_comm_grid,ier)
           bufind=0
           do jatm=1,natom
             jtypat=typat(atindx1(jatm))
             lm_sizej=pawtab(jtypat)%lcut_size**2
             if (ngrad> 0) prod_nondiag(jatm)%value(:,:)= &
&             buf(1:ngrad_nondiag,bufind+1:bufind+lm_sizej)
             if (ngradp>0) prodp_nondiag(jatm)%value(:,:)= &
&             buf(ngrad_nondiag+1:ngrad_nondiag+ngradp_nondiag,bufind+1:bufind+lm_sizej)
             bufind=bufind+lm_sizej*(ngrad_nondiag+ngradp_nondiag)
           end do
           ABI_DEALLOCATE(buf)
         end if
       end if

!      ----------------------------------------------------------------
!      Compute final sums (i.e. derivatives of Sum_ij[rho_ij.Intg{Qij.Vloc}]

!      ---- Compute terms common to all gradients
       jrhoij=1
       do irhoij=1,pawrhoij_iatom%nrhoijsel
         klmn=pawrhoij_iatom%rhoijselect(irhoij)
         klm =pawtab(itypat)%indklmn(1,klmn)
         lmin=pawtab(itypat)%indklmn(3,klmn)
         lmax=pawtab(itypat)%indklmn(4,klmn)
         ro =pawrhoij_iatom%rhoijp(jrhoij,ispden)
         ro_d=ro*pawtab(itypat)%dltij(klmn)
         do ll=lmin,lmax,2
           do ilm=ll**2+1,(ll+1)**2
             isel=pawang%gntselect(ilm,klm)
             if (isel>0) then
               grhat_x=ro_d*pawtab(itypat)%qijl(ilm,klmn)
               do mu=1,ngrad
                 grhat_tmp(mu,idiag)=grhat_tmp(mu,idiag)+grhat_x*prod(mu,ilm)
               end do
             end if
           end do
         end do
         jrhoij=jrhoij+pawrhoij_iatom%cplex_rhoij
       end do

!      ---- Add additional (diagonal) terms for dynamical matrix
!      ---- Terms including rhoij derivatives
       if (optgr2==1) then
         klmn1=1
         do klmn=1,lmn2_size
           klm =pawtab(itypat)%indklmn(1,klmn)
           lmin=pawtab(itypat)%indklmn(3,klmn)
           lmax=pawtab(itypat)%indklmn(4,klmn)
           dlt_tmp=pawtab(itypat)%dltij(klmn)
           do ll=lmin,lmax,2
             do ilm=ll**2+1,(ll+1)**2
               isel=pawang%gntselect(ilm,klm)
               if (isel>0) then
                 ro_d= dlt_tmp*pawtab(itypat)%qijl(ilm,klmn)
                 do mu=1,9
                   mua=alpha(mu);mub=beta(mu)
                   grhat_tmp(ishift_gr2+mu,idiag)=grhat_tmp(ishift_gr2+mu,idiag)&
&                   +ro_d*pawrhoij_iatom%grhoij(ishift_grhoij+mua,klmn1,ispden)*prodp(mub+ishift_gr,ilm)
                 end do
               end if
             end do
           end do
           klmn1=klmn1+pawrhoij_iatom%cplex_rhoij
         end do ! klmn
       end if ! optgr2

!      ---- Add additional (diagonal) terms for elastic tensor
!      ---- Terms including rhoij derivatives
       if (optstr2==1)then
         klmn1=1
         do klmn=1,lmn2_size
           klm =pawtab(itypat)%indklmn(1,klmn)
           lmin=pawtab(itypat)%indklmn(3,klmn)
           lmax=pawtab(itypat)%indklmn(4,klmn)
           dlt_tmp=pawtab(itypat)%dltij(klmn)
           do ll=lmin,lmax,2
             do ilm=ll**2+1,(ll+1)**2
               isel=pawang%gntselect(ilm,klm)
               if (isel>0) then
                 ro_d=dlt_tmp*pawtab(itypat)%qijl(ilm,klmn)
                 mu=1
                 do mua=1,6
                   do mub=1,6
                     grhat_tmp(ishift_str2+mu,iatm)= grhat_tmp(ishift_str2+mu,iatm)&
&                     +ro_d*pawrhoij_iatom%grhoij(mub,klmn1,ispden)*prodp(mua,ilm)
                     grhat_tmp(ishift_str2+mu,iatm)= grhat_tmp(ishift_str2+mu,iatm)&
&                     +ro_d*pawrhoij_iatom%grhoij(mua,klmn1,ispden)*prodp(mub,ilm)
                     mu=mu+1
                   end do
!                  INTERNAL STRAIN CONTRIBUTION
                   do idir=1,3
                     grhat_tmp(ishift_str2is+(mua-1)*3+idir,iatm) = grhat_tmp(ishift_str2is+(mua-1)*3+idir,iatm)&
&                     +ro_d*pawrhoij_iatom%grhoij(ishift_grhoij+idir,klmn1,ispden)*prodp(mua,ilm)
                     grhat_tmp(ishift_str2is+(mua-1)*3+idir,iatm) = grhat_tmp(ishift_str2is+(mua-1)*3+idir,iatm)&
&                     +ro_d*pawrhoij_iatom%grhoij(mua,klmn1,ispden)*prodp(6+idir,ilm)
                   end do
                 end do
               end if
             end do
           end do
           klmn1=klmn1+pawrhoij_iatom%cplex_rhoij
         end do
       end if ! optstr2

!      ---- Add off-diagonal additional contributions for second gradients
       if (optgr2==1.or.optstr2==1) then
         do jatm=1,natom
           jatom_tot=atindx1(jatm);jtypat=typat(jatom_tot)
           pawrhoij_jatom => pawrhoij_tot(jatom_tot)

!          ---- Dynamical matrix
           if (optgr2==1) then

!            Off-diagonal term including rhoij
             if (dyfr_cplex==1.or.cplex==1) then
               jrhoij=1
               do irhoij=1,pawrhoij_jatom%nrhoijsel
                 klmn=pawrhoij_jatom%rhoijselect(irhoij)
                 klm =pawtab(jtypat)%indklmn(1,klmn)
                 lmin=pawtab(jtypat)%indklmn(3,klmn)
                 lmax=pawtab(jtypat)%indklmn(4,klmn)
                 ro  =pawrhoij_jatom%rhoijp(jrhoij,ispden)
                 ro_d=ro*pawtab(jtypat)%dltij(klmn)
                 do ll=lmin,lmax,2
                   do ilm=ll**2+1,(ll+1)**2
                     isel=pawang%gntselect(ilm,klm)
                     if (isel>0) then
                       grhat_x=ro_d*pawtab(jtypat)%qijl(ilm,klmn)
                       do mu=1,9
                         grhat_tmp(ishift_gr2+mu,jatm)=grhat_tmp(ishift_gr2+mu,jatm) &
&                         +grhat_x*prod_nondiag(jatm)%value(ishift_gr2+mu,ilm)
                       end do
                     end if
                   end do
                 end do
                 jrhoij=jrhoij+pawrhoij_jatom%cplex_rhoij
               end do
             else
               jrhoij=1;mushift=ishift_gr2+9
               do irhoij=1,pawrhoij_jatom%nrhoijsel
                 klmn=pawrhoij_jatom%rhoijselect(irhoij)
                 klm =pawtab(jtypat)%indklmn(1,klmn)
                 lmin=pawtab(jtypat)%indklmn(3,klmn)
                 lmax=pawtab(jtypat)%indklmn(4,klmn)
                 ro  =pawrhoij_jatom%rhoijp(jrhoij,ispden)
                 ro_d=ro*pawtab(jtypat)%dltij(klmn)
                 do ll=lmin,lmax,2
                   do ilm=ll**2+1,(ll+1)**2
                     isel=pawang%gntselect(ilm,klm)
                     if (isel>0) then
                       grhat_x=ro_d*pawtab(jtypat)%qijl(ilm,klmn)
                       do mu=1,9
                         grhat_tmp(ishift_gr2+mu,jatm)=grhat_tmp(ishift_gr2+mu,jatm)&
&                         +grhat_x*prod_nondiag(jatm)%value(ishift_gr2+mu,ilm)
                         grhat_tmp(mushift+mu,jatm)=grhat_tmp(mushift+mu,jatm)&
&                         +grhat_x*prod_nondiag(jatm)%value(ishift_gr2+9+mu,ilm)
                       end do
                     end if
                   end do
                 end do
                 jrhoij=jrhoij+pawrhoij_jatom%cplex_rhoij
               end do
             end if

!            Off-diagonal term including rhoij derivative
             if (dyfr_cplex==1.or.cplex==1) then
               klmn1=1
               do klmn=1,pawrhoij_jatom%lmn2_size
                 klm =pawtab(jtypat)%indklmn(1,klmn)
                 lmin=pawtab(jtypat)%indklmn(3,klmn)
                 lmax=pawtab(jtypat)%indklmn(4,klmn)
                 dlt_tmp=pawtab(jtypat)%dltij(klmn)
                 do ll=lmin,lmax,2
                   do ilm=ll**2+1,(ll+1)**2
                     isel=pawang%gntselect(ilm,klm)
                     if (isel>0) then
                       ro_d=dlt_tmp*pawtab(jtypat)%qijl(ilm,klmn)
                       do mu=1,9
                         mua=alpha(mu);mub=beta(mu)
                         grhat_tmp(ishift_gr2+mu,jatm)=grhat_tmp(ishift_gr2+mu,jatm) &
&                         +ro_d*pawrhoij_jatom%grhoij(ishift_grhoij+mua,klmn1,ispden) &
&                         *prodp_nondiag(jatm)%value(ishift2_gr+mub,ilm)
                       end do
                     end if
                   end do
                 end do
                 klmn1=klmn1+pawrhoij_jatom%cplex_rhoij
               end do ! klmn
             else ! ngradp_nondiag>=6
               klmn1=1;mushift=ishift_gr2+9
               do klmn=1,pawrhoij_jatom%lmn2_size
                 klm =pawtab(jtypat)%indklmn(1,klmn)
                 lmin=pawtab(jtypat)%indklmn(3,klmn)
                 lmax=pawtab(jtypat)%indklmn(4,klmn)
                 dlt_tmp=pawtab(jtypat)%dltij(klmn)
                 do ll=lmin,lmax,2
                   do ilm=ll**2+1,(ll+1)**2
                     isel=pawang%gntselect(ilm,klm)
                     if (isel>0) then
                       ro_d=dlt_tmp*pawtab(jtypat)%qijl(ilm,klmn)
                       do mu=1,9
                         mua=alpha(mu);mub=beta(mu)
                         grhat_tmp(ishift_gr2+mu,jatm)=grhat_tmp(ishift_gr2+mu,jatm) &
&                         +ro_d*pawrhoij_jatom%grhoij(ishift_grhoij+mua,klmn1,ispden) &
&                         *prodp_nondiag(jatm)%value(ishift2_gr+mub,ilm)
                         grhat_tmp(mushift+mu,jatm)=grhat_tmp(mushift+mu,jatm) &
&                         +ro_d*pawrhoij_jatom%grhoij(ishift_grhoij+mua,klmn1,ispden) &
&                         *prodp_nondiag(jatm)%value(ishift2_gr+3+mub,ilm)
                       end do
                     end if
                   end do
                 end do
                 klmn1=klmn1+pawrhoij_jatom%cplex_rhoij
               end do
             end if
           end if ! optgr2

!          ---- Elastic tensor
           if (optstr2==1)then

!            Off-diagonal term including rhoij
             jrhoij=1;
             do irhoij=1,pawrhoij_jatom%nrhoijsel
               klmn=pawrhoij_jatom%rhoijselect(irhoij)
               klm =pawtab(jtypat)%indklmn(1,klmn)
               lmin=pawtab(jtypat)%indklmn(3,klmn)
               lmax=pawtab(jtypat)%indklmn(4,klmn)
               ro  =pawrhoij_jatom%rhoijp(jrhoij,ispden)
               ro_d=ro*pawtab(jtypat)%dltij(klmn)
               do ll=lmin,lmax,2
                 do ilm=ll**2+1,(ll+1)**2
                   isel=pawang%gntselect(ilm,klm)
                   if (isel>0) then
                     grhat_x=ro_d*pawtab(jtypat)%qijl(ilm,klmn)
                     do mu=1,18
                       grhat_tmp2(mu,jatm)=grhat_tmp2(mu,jatm) &
&                       +grhat_x*prod_nondiag(jatm)%value(ishift_str2is+mu,ilm)
                     end do
                   end if
                 end do
               end do
               jrhoij=jrhoij+pawrhoij_jatom%cplex_rhoij
             end do
!            Off-diagonal term including rhoij derivative
             klmn1=1
             do klmn=1,pawrhoij_jatom%lmn2_size
               klm =pawtab(jtypat)%indklmn(1,klmn)
               lmin=pawtab(jtypat)%indklmn(3,klmn)
               lmax=pawtab(jtypat)%indklmn(4,klmn)
               dlt_tmp=pawtab(jtypat)%dltij(klmn)
               do ll=lmin,lmax,2
                 do ilm=ll**2+1,(ll+1)**2
                   isel=pawang%gntselect(ilm,klm)
                   if (isel>0) then
                     ro_d=dlt_tmp*pawtab(jtypat)%qijl(ilm,klmn)
                     mu=1
                     do mua=1,6
                       do idir=1,3
                         grhat_tmp2((mua-1)*3+idir,jatm) = grhat_tmp2((mua-1)*3+idir,jatm) &
&                         +ro_d*pawrhoij_jatom%grhoij(mua,klmn1,ispden) &
&                         *prodp_nondiag(jatm)%value(idir,ilm)
                       end do
                     end do
                   end if
                 end do
               end do
               klmn1=klmn1+pawrhoij_jatom%cplex_rhoij
             end do
           end if ! optstr2

         end do ! jatm
       end if ! optgr2 or optstr2

!    ----------------------------------------------------------------
!    End of loop over spin components

     end do ! ispden

!    Eventually free temporary space for g_l(r).Y_lm(r) factors
     if (pawfgrtab_iatom%gylm_allocated==2) then
       ABI_DEALLOCATE(pawfgrtab_iatom%gylm)
       ABI_ALLOCATE(pawfgrtab_iatom%gylm,(0,0))
       pawfgrtab_iatom%gylm_allocated=0
     end if
     if (pawfgrtab_iatom%gylmgr_allocated==2) then
       ABI_DEALLOCATE(pawfgrtab_iatom%gylmgr)
       ABI_ALLOCATE(pawfgrtab_iatom%gylmgr,(0,0,0))
       pawfgrtab_iatom%gylmgr_allocated=0
     end if
     if (pawfgrtab_iatom%gylmgr2_allocated==2) then
       ABI_DEALLOCATE(pawfgrtab_iatom%gylmgr2)
       ABI_ALLOCATE(pawfgrtab_iatom%gylmgr2,(0,0,0))
       pawfgrtab_iatom%gylmgr2_allocated=0
     end if
     if (pawfgrtab_iatom%expiqr_allocated==2) then
       ABI_DEALLOCATE(pawfgrtab_iatom%expiqr)
       ABI_ALLOCATE(pawfgrtab_iatom%expiqr,(0,0))
       pawfgrtab_iatom%expiqr_allocated=0
     end if

!    ----------------------------------------------------------------
!    Copy results in corresponding arrays

!    ==== Forces ====
!    Convert from cartesian to reduced coordinates
     if (optgr==1) then
       mushift=3*(iatm-1)
       tmp(1:3)=grhat_tmp(ishift_gr+1:ishift_gr+3,idiag)
       do mu=1,3
         hatgr(mu+mushift)=rprimd(1,mu)*tmp(1)+rprimd(2,mu)*tmp(2)+rprimd(3,mu)*tmp(3)
       end do
     end if

!    ==== Stresses ====
     if (optstr==1) then
       hatstr(1:6)=hatstr(1:6)+grhat_tmp(ishift_str+1:ishift_str+6,idiag)
     end if

!    ==== Frozen wf part of dyn. matrix ====
     if (optgr2==1) then
       do jatm=1,natom
         do mu=1,9
           mua=alpha(mu);mub=beta(mu)
           dyfr(1,mub,mua,jatm,iatm)=grhat_tmp(ishift_gr2+mu,jatm)
         end do
         if (dyfr_cplex==2.and.cplex==2) then
           mushift=ishift_gr2+9
           do mu=1,9
             mua=alpha(mu);mub=beta(mu)
             dyfr(2,mub,mua,jatm,iatm)=grhat_tmp(mushift+mu,jatm)
           end do
         end if
       end do
     end if

!    ==== Elastic tensor ====
     if (optstr2==1) then
       eltfr(1:6,1:6)=eltfr(1:6,1:6)+ &
&       reshape(grhat_tmp(ishift_str2+1:ishift_str2+36,iatm),(/6,6/))
!      Convert internal Strain in reduced coordinates
       do mua = 1,6
         tmp(1:3)=grhat_tmp(ishift_str2is+(mua-1)*3+1:ishift_str2is+(mua-1)*3+3,iatm)
         do idir=1,3
           eltfr(6+(iatm-1)*3+idir,mua)=eltfr(6+(iatm-1)*3+idir,mua)+ &
&           (rprimd(1,idir)*tmp(1)+rprimd(2,idir)*tmp(2)+rprimd(3,idir)*tmp(3))
         end do
         do jatm=1,natom
           tmp(1:3)=grhat_tmp2((mua-1)*3+1:(mua-1)*3+3,jatm)
           do idir=1,3
             eltfr(6+(iatm-1)*3+idir,mua)=eltfr(6+(iatm-1)*3+idir,mua)+ &
&             (rprimd(1,idir)*tmp(1)+rprimd(2,idir)*tmp(2)+rprimd(3,idir)*tmp(3))
           end do
         end do
       end do
     end if

!    ----------------------------------------------------------------
!    End loops on types and atoms

     ABI_DEALLOCATE(vloc)
     if (ngrad>0)  then
       ABI_DEALLOCATE(prod)
     end if
     if (ngradp>0)  then
       ABI_DEALLOCATE(prodp)
     end if
     if (optgr2==1.or.optstr2==1) then
       do jatm=1,natom
         ABI_DEALLOCATE(prod_nondiag(jatm)%value)
         ABI_DEALLOCATE(prodp_nondiag(jatm)%value)
       end do
     end if
   end do
   iatshft=iatshft+nattyp(itypat)
 end do

!Reduction in case of parallelisation over atoms
 if (paral_atom) then
   bufsiz=3*natom*optgr+6*optstr
   if (save_memory) bufsiz=bufsiz+9*dyfr_cplex*natom**2*optgr2+6*(6+3*natom)*optstr2
   if (bufsiz>0) then
     ABI_ALLOCATE(buf1,(bufsiz))
     if (optgr==1) buf1(1:3*natom)=hatgr(1:3*natom)
     indx=optgr*3*natom
     if (optstr==1) buf1(indx+1:indx+6)=hatstr(1:6)
     indx=indx+optstr*6
     if (save_memory) then
       if (optgr2==1) then
         buf1(indx+1:indx+9*dyfr_cplex*natom**2)= &
&         reshape(dyfr,(/9*dyfr_cplex*natom**2/))
         indx=indx+9*dyfr_cplex*natom**2
       end if
       if (optstr2==1) then
         buf1(indx+1:indx+6*(6+3*natom))= &
&         reshape(eltfr,(/6*(6+3*natom)/))
         indx=indx+6*(6+3*natom)
       end if
     end if
     call xmpi_sum(buf1,my_comm_atom,ier)
     if (optgr==1) hatgr(1:3*natom)=buf1(1:3*natom)
     indx=optgr*3*natom
     if (optstr==1) hatstr(1:6)=buf1(indx+1:indx+6)
     indx=indx+optstr*6
     if (save_memory) then
       if (optgr2==1) then
         dyfr(1:dyfr_cplex,1:3,1:3,1:natom,1:natom)= &
&         reshape(buf1(indx+1:indx+9*dyfr_cplex*natom**2),(/dyfr_cplex,3,3,natom,natom/))
         indx=indx+9*dyfr_cplex*natom**2
       end if
       if (optstr2==1) then
         eltfr(1:6+3*natom,1:6)= &
&         reshape(buf1(indx+1:indx+6*(6+3*natom)),(/6+3*natom,6/))
         indx=indx+6*(6+3*natom)
       end if
     end if
     ABI_DEALLOCATE(buf1)
   end if
 end if

!Deallocate additional memory
 ABI_DEALLOCATE(grhat_tmp)
 if (optgr2==1.or.optstr2==1) then
   ABI_DEALLOCATE(mu4)
   ABI_DEALLOCATE(atindx)
   if (optgr2==1.or.optstr2==1) then
     ABI_DEALLOCATE(vpsp1_gr)
   end if
   if (optstr2==1) then
     ABI_DEALLOCATE(grhat_tmp2)
     ABI_DEALLOCATE(vpsp1_str)
   end if
   ABI_DATATYPE_DEALLOCATE(prod_nondiag)
   ABI_DATATYPE_DEALLOCATE(prodp_nondiag)
   if (.not.save_memory) then
     do jatom=1,size(pawfgrtab)
       pawfgrtab_jatom => pawfgrtab(jatom)
       if (pawfgrtab(jatom)%gylm_allocated==2) then
         ABI_DEALLOCATE(pawfgrtab(jatom)%gylm)
         ABI_ALLOCATE(pawfgrtab(jatom)%gylm,(0,0))
         pawfgrtab(jatom)%gylm_allocated=0
       end if
       if (pawfgrtab(jatom)%gylmgr_allocated==2) then
         ABI_DEALLOCATE(pawfgrtab(jatom)%gylmgr)
         ABI_ALLOCATE(pawfgrtab(jatom)%gylmgr,(0,0,0))
         pawfgrtab(jatom)%gylmgr_allocated=0
       end if
       if (pawfgrtab(jatom)%gylmgr2_allocated==2) then
         ABI_DEALLOCATE(pawfgrtab(jatom)%gylmgr2)
         ABI_ALLOCATE(pawfgrtab(jatom)%gylmgr2,(0,0,0))
         pawfgrtab(jatom)%gylmgr2_allocated=0
       end if
       if (pawfgrtab(jatom)%expiqr_allocated==2) then
         ABI_DEALLOCATE(pawfgrtab(jatom)%expiqr)
         ABI_ALLOCATE(pawfgrtab(jatom)%expiqr,(0,0))
         pawfgrtab(jatom)%expiqr_allocated=0
       end if
     end do
   end if
   if (paral_atom) then
     if ((.not.save_memory).and.paral_atom_pawfgrtab) then
       call pawfgrtab_free(pawfgrtab_tot)
       ABI_DATATYPE_DEALLOCATE(pawfgrtab_tot)
     end if
     if (paral_atom_pawrhoij) then
       call pawrhoij_free(pawrhoij_tot)
       ABI_DATATYPE_DEALLOCATE(pawrhoij_tot)
     end if
   end if
 end if

!----------------------------------------------------------------------
!Update non-local gardients

!===== Update forces =====
 if (optgr==1) then
   grnl(1:3*natom)=grnl(1:3*natom)+hatgr(1:3*natom)
   ABI_DEALLOCATE(hatgr)
 end if

!===== Convert stresses (add diag and off-diag contributions) =====
 if (optstr==1) then

!  Has to compute int[nhat*vtrial]
   hatstr_diag=zero
   if (nspden==1.or.dimvtrial==1) then
     do ic=1,nfft
       hatstr_diag=hatstr_diag+vtrial_(ic,1)*nhat(ic,1)
     end do
   else if (nspden==2) then
     do ic=1,nfft
       hatstr_diag=hatstr_diag+vtrial_(ic,1)*nhat(ic,2)+vtrial_(ic,2)*(nhat(ic,1)-nhat(ic,2))
     end do
   else if (nspden==4) then
     do ic=1,nfft
       hatstr_diag=hatstr_diag+half*(vtrial_(ic,1)*(nhat(ic,1)+nhat(ic,4)) &
&       +vtrial_(ic,2)*(nhat(ic,1)-nhat(ic,4))) &
&       +vtrial_(ic,3)*nhat(ic,2)+vtrial_(ic,4)*nhat(ic,3)
     end do
   end if
   hatstr_diag=hatstr_diag*fact_ucvol
   if (paral_grid) then
     call xmpi_sum(hatstr_diag,my_comm_grid,ier)
   end if

!  Convert hat contribution

   hatstr(1:3)=(hatstr(1:3)+hatstr_diag)/ucvol
   hatstr(4:6)= hatstr(4:6)/ucvol

!  Add to already computed NL contrib
   nlstr(1:6)=nlstr(1:6)+hatstr(1:6)

!  Apply symmetries
   call stresssym(gprimd,nsym,nlstr,symrec)
 end if

!===== Convert dynamical matrix (from cartesian to reduced coordinates) =====
 if (optgr2==1) then
   do iatm=1,natom
     do jatm=1,natom
       do mua=1,3
         do mub=1,3
           work1(1,mua,mub)=dyfr(1,mub,mua,jatm,iatm)+dyfr(1,mua,mub,iatm,jatm)
         end do
       end do
       if (dyfr_cplex==2) then
         do mua=1,3
           do mub=1,3
             work1(2,mua,mub)=dyfr(2,mub,mua,jatm,iatm)-dyfr(2,mua,mub,iatm,jatm)
           end do
         end do
       end if
       do mu=1,3
         work2(:,:,mu)=rprimd(1,mu)*work1(:,:,1)+rprimd(2,mu)*work1(:,:,2)+rprimd(3,mu)*work1(:,:,3)
       end do
       do mub=1,3
         do mua=1,3
           dyfrnl(:,mua,mub,jatm,iatm)=dyfrnl(:,mua,mub,jatm,iatm) &   ! Already contains NL projectors contribution
&          +rprimd(1,mua)*work2(:,1,mub) &
&           +rprimd(2,mua)*work2(:,2,mub) &
&           +rprimd(3,mua)*work2(:,3,mub)
         end do
       end do
     end do
   end do
   ABI_DEALLOCATE(dyfr)
 end if

!===== Update elastic tensor =====
 if (optstr2==1) then
   eltfrnl(1:6+3*natom,1:6)=eltfrnl(1:6+3*natom,1:6)+eltfr(1:6+3*natom,1:6)
   ABI_DEALLOCATE(eltfr)
 end if

!----------------------------------------------------------------------
!End

!Destroy temporary space
 if (usexcnhat==0)  then
   ABI_DEALLOCATE(vtrial_)
 end if

!Destroy atom tables used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)
 if (paral_atom) then
   ABI_DEALLOCATE(atm_indx)
 end if

!Destroy FFT tables used for parallelism
 if ((optgr2==1.or.optstr2==1).and.(.not.present(comm_fft))) then
   call destroy_distribfft(my_distribfft)
   ABI_DATATYPE_DEALLOCATE(my_distribfft)
 end if

 DBG_ENTER("COLL")

 CONTAINS
!!***

! ------------------------------------------------
!!****f* pawgrnl/pawgrnl_convert
!! NAME
!!  pawgrnl_convert
!!
!! FUNCTION
!!  notation: Convert index of the elastic tensor:
!!    - voigt notation       => 32
!!    - normal notation      => 3 3 2 2
!!    - notation for gylmgr2 => 32 32 32 32 => 4 4 4
!!
!! INPUTS
!!  eps_alpha, eps_beta, eps_delta, eps_gamma
!!
!! OUTPUT
!!  mu4(4) = array with index for the second derivative of gylm
!!
!! SIDE EFFECTS
!!  mu4(4) = input : array with index for the second derivative of gylm
!!           output: the 4 indexes for the calculation of the second derivative of gylm
!!
!! PARENTS
!!      m_paw_dfpt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawgrnl_convert(mu4,eps_alpha,eps_beta,eps_gamma,eps_delta)

!Arguments ------------------------------------
 !scalar
 integer,intent(in)  :: eps_alpha,eps_beta
 integer,optional,intent(in)  :: eps_gamma,eps_delta
 !array
 integer,intent(inout) :: mu4(4)

!Local variables-------------------------------
 integer :: eps1,eps2,i,j,k
 integer,allocatable :: mu_temp(:)

! *************************************************************************

 ABI_ALLOCATE(mu_temp,(4))
 if (present(eps_gamma).and.present(eps_delta)) then
   mu_temp(1)=eps_alpha
   mu_temp(2)=eps_beta
   mu_temp(3)=eps_gamma
   mu_temp(4)=eps_delta
 else
   mu_temp(1)=eps_alpha
   mu_temp(2)=eps_beta
   mu_temp(3)= 0
   mu_temp(4)= 0
 end if
 k=1
 do i=1,2
   eps1=mu_temp(i)
   do j=1,2
     eps2=mu_temp(2+j)
     if(eps1==eps2) then
       if(eps1==1) mu4(k)=1;
       if(eps1==2) mu4(k)=2;
       if(eps1==3) mu4(k)=3;
     else
       if((eps1==3.and.eps2==2).or.(eps1==2.and.eps2==3)) mu4(k)=4;
       if((eps1==3.and.eps2==1).or.(eps1==1.and.eps2==3)) mu4(k)=5;
       if((eps1==1.and.eps2==2).or.(eps1==2.and.eps2==1)) mu4(k)=6;
     end if
     k=k+1
   end do
 end do
 ABI_DEALLOCATE(mu_temp)

end subroutine pawgrnl_convert
! ------------------------------------------------

end subroutine pawgrnl
!!***

!----------------------------------------------------------------------

!!****f* m_paw_dfpt/dsdr_k_paw
!! NAME
!! dsdr_k_paw
!!
!! FUNCTION
!! compute on-site terms for forces and stresses for finite electric fields with PAW
!!
!! INPUTS
!!  cprj_k (pawcprj_type) :: cprj for occupied bands at point k
!!  cprj_kb :: cprj for occupied bands at point k+b
!!  dtefield :: structure referring to all efield and berry's phase variables
!!  kdir :: integer giving direction along which overlap is computed for ket
!!  kfor :: integer indicating whether to compute forward (1) or backward (2)
!!    along kpt string
!!  natom :: number of atoms in cell
!!  typat :: typat(natom) type of each atom
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! dsdr :: array of the on-site PAW parts of the derivatives with respect to atm
!!               positions and/or strains of the overlaps between Bloch states at points
!!               k and k+b, for the various pairs of bands
!!
!! NOTES
!! This routine assumes that the cprj are not explicitly ordered by
!! atom type.
!!
!! PARENTS
!!      m_berryphase_new
!!
!! CHILDREN
!!
!! SOURCE

 subroutine dsdr_k_paw(cprj_k,cprj_kb,dsdr,dtefield,kdir,kfor,mband,natom,ncpgr,typat)

!Arguments---------------------------
!scalars
 integer,intent(in) :: kdir,kfor,mband,natom,ncpgr
 character(len=500) :: message
 type(efield_type),intent(in) :: dtefield
 type(pawcprj_type),intent(in) :: cprj_k(natom,dtefield%nspinor*mband)
 type(pawcprj_type),intent(in) :: cprj_kb(natom,dtefield%nspinor*mband)

!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(inout) :: dsdr(2,natom,ncpgr,dtefield%mband_occ,dtefield%mband_occ)

!Local variables---------------------------
!scalars
 integer :: iatom,iband,ibs,icpgr,ilmn,ispinor,itypat
 integer :: jband,jbs,jlmn,klmn,nspinor
 complex(dpc) :: cpk,cpkb,dcpk,dcpkb,cterm,paw_onsite

! *************************************************************************

!initialize dsdr
 dsdr(:,:,:,:,:) = zero

! if 3 gradients we are in the ctocprj choice 2 case
! and the 3 gradients are due to the atomic displacements
! if 6 gradients we are in the ctocprj choice 3 case
! and the 6 gradients are due to the strains
! if 9 gradients we are in the ctocprj choice 23 case
! and the first six are due to strain, last three due to displacements
 if (ncpgr /= 3 .and. ncpgr /= 6 .and. ncpgr /= 9) then
   message = ' dsdr_k_paw called with ncpgr /= 3, 6, or 9 (no gradients) '
   ABI_BUG(message)
 end if

 nspinor = dtefield%nspinor

 do iatom = 1, natom
   itypat = typat(iatom)

   do ilmn=1,dtefield%lmn_size(itypat)
     do jlmn=1,dtefield%lmn_size(itypat)
       klmn=max(ilmn,jlmn)*(max(ilmn,jlmn)-1)/2 + min(ilmn,jlmn)
       paw_onsite = cmplx(dtefield%qijb_kk(1,klmn,iatom,kdir),&
&       dtefield%qijb_kk(2,klmn,iatom,kdir))
       if (kfor > 1) paw_onsite = conjg(paw_onsite)
       do iband = 1, dtefield%mband_occ
         do jband = 1, dtefield%mband_occ
           do ispinor = 1, nspinor
             do icpgr = 1, ncpgr
               ibs = nspinor*(iband-1) + ispinor
               jbs = nspinor*(jband-1) + ispinor
               cpk=cmplx(cprj_k(iatom,ibs)%cp(1,ilmn),cprj_k(iatom,ibs)%cp(2,ilmn))
               dcpk=cmplx(cprj_k(iatom,ibs)%dcp(1,icpgr,ilmn),cprj_k(iatom,ibs)%dcp(2,icpgr,ilmn))
               cpkb=cmplx(cprj_kb(iatom,jbs)%cp(1,jlmn),cprj_kb(iatom,jbs)%cp(2,jlmn))
               dcpkb=cmplx(cprj_kb(iatom,jbs)%dcp(1,icpgr,jlmn),cprj_kb(iatom,jbs)%dcp(2,icpgr,jlmn))
               cterm=paw_onsite*(conjg(dcpk)*cpkb+conjg(cpk)*dcpkb)
               dsdr(1,iatom,icpgr,iband,jband) = dsdr(1,iatom,icpgr,iband,jband)+real(cterm)
               dsdr(2,iatom,icpgr,iband,jband) = dsdr(2,iatom,icpgr,iband,jband)+aimag(cterm)
             end do ! end loop over icpgr
           end do ! end loop over ispinor
         end do ! end loop over jband
       end do ! end loop over iband
     end do ! end loop over ilmn
   end do ! end loop over jlmn

 end do ! end loop over atoms

 end subroutine    dsdr_k_paw
!!***

!----------------------------------------------------------------------

END MODULE m_paw_dfpt
!!***

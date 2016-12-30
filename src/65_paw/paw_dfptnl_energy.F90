!{\src2tex{textfont=tt}}
!!****f* ABINIT/paw_dfptnl_energy
!! NAME
!! paw_dfptnl_energy
!!
!! FUNCTION
!! This routine compute the Hartree+XC PAW on-site contributions to a 1st-order or 2nd-order energy.
!!  These contributions are equal to:
!!    E_onsite=
!!       Int{ VHxc[n1_a^(1);nc^(1)].n1_b }
!!      -Int{ VHxc[tild_n1_a^(j1)+hat_n1_a^(j1);tild_n_c^(j1)].(tild_n1_b+n1_b)^(j2) }
!! Some typical uses:
!!  A-Contribution to non-stationary expression of the 2nd-order total energy:
!!    In that case, n1_a^(1)[r]=n1^(j1)[r] and n1_b[r]=delta_n1^(j2)[r]
!!    where j1 and j2 are two given perturbations,
!!    and delta_n1^(j)[r] is the 1s-order density only due to change of WF overlap.
!!    See PRB 78, 035105 (2008), Eq.(80)
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
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!!
!! PARENTS
!!      dfpt_nstpaw,newfermie1
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,pawdensities,pawdijhartree,pawxc3,pawxcm3
!!      timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine paw_dfptnl_energy(d3exc,ixc,my_natom,natom,ntypat,&
&                    paw_an0,pawang,pawprtvol,pawrad,&
&                    pawrhoij_1,pawrhoij_2,pawrhoij_3,&
&                    pawtab,pawxcdev,&
&                    mpi_atmtab,comm_atom) ! optional arguments (parallelism)


 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi, only : xmpi_comm_self,xmpi_sum

 use m_pawang,     only : pawang_type
 use m_pawrad,     only : pawrad_type
 use m_pawtab,     only : pawtab_type
 use m_paw_an,     only : paw_an_type
 use m_paw_ij,     only : paw_ij_type
 use m_pawrhoij,   only : pawrhoij_type
 use m_pawdij,     only : pawdijhartree
 use m_pawxc,      only : pawxc3, pawxcm3
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_dfptnl_energy'
 use interfaces_65_paw, except_this_one => paw_dfptnl_energy
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ixc,my_natom,natom,ntypat
 integer,intent(in) :: pawprtvol,pawxcdev!,xclevel
 integer,optional,intent(in) :: comm_atom
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(out) :: d3exc(2)
 type(paw_an_type),intent(in) :: paw_an0(my_natom)
! type(paw_an_type),intent(inout) :: paw_an1(my_natom)
! type(paw_ij_type),intent(inout) :: paw_ij1(my_natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij_1(natom)
 type(pawrhoij_type),intent(in) :: pawrhoij_2(natom)
 type(pawrhoij_type),intent(in) :: pawrhoij_3(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: cplex_1,cplex_2,cplex_3,cplex_dijh1,iatom,iatom_tot,ierr,ipts,irhoij,ispden,itypat,jrhoij
 integer :: klmn,ilm,lm_size_all,lm_size_eff,mesh_size,my_comm_atom,npts,nspden,nspdiag,nzlmopt
 integer :: opt_compch,optexc,optvxc,usecore,usetcore,usexcnhat
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: compch,d3exc1_iat(2),eexc,eexc_im
 character(len=500) :: msg
!arrays
 integer,pointer :: my_atmtab(:)
 logical,allocatable :: lmselect_1(:),lmselect_2(:),lmselect_3(:),lmselect_tmp(:)
 real(dp) :: dij(2),delta_energy_h(2),delta_energy_xc(2),ro(2),tsec(2)
 real(dp),allocatable :: nhat1_1(:,:,:),rho1_1(:,:,:),trho1_1(:,:,:) ! kxc_dum(:,:,:),
 real(dp),allocatable :: nhat1_2(:,:,:),rho1_2(:,:,:),trho1_2(:,:,:)
 real(dp),allocatable :: nhat1_3(:,:,:),rho1_3(:,:,:),trho1_3(:,:,:)
 real(dp),allocatable :: rho1arr(:,:),rho2arr(:,:),rho3arr(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 nzlmopt = 0

 if (pawxcdev/=0) then
   MSG_BUG("paw_dfptnl_energy is not implemented for pawxcdev /=0")
 end if

! call timab(567,1,tsec)

! if (.not.(ipert1==natom+1.or.ipert1==natom+10.or.ipert1==natom+11 &
!& .or.ipert2==natom+1.or.ipert2==natom+10.or.ipert2==natom+11)) then
!   if((abs(nzlmopt_a)/=1.and.nzlmopt_a/=0).or.(abs(nzlmopt_b)/=1.and.nzlmopt_b/=0)) then
!     msg='invalid value for nzlmopt !'
!     MSG_BUG(msg)
!   end if
!   if (my_natom>0) then
!     if(paw_ij1(1)%has_dijhartree==0) then
!       msg='dijhartree must be allocated !'
!       MSG_BUG(msg)
!     end if
!     if(paw_an1(1)%has_vxc==0) then
!       msg='vxc1 and vxct1 must be allocated !'
!       MSG_BUG(msg)
!     end if
!     if(paw_an0(1)%has_kxc==0) then
!       msg='kxc1 must be allocated !'
!       MSG_BUG(msg)
!     end if
!     if ((ipert1<=natom.or.ipert1==natom+1.or.ipert1==natom+10.or.ipert1==natom+11).and.paw_an0(1)%has_kxc/=2) then
!       msg='XC kernels for ground state must be in memory !'
!       MSG_BUG(msg)
!     end if
!     if (paw_ij1(1)%cplex/=paw_an1(1)%cplex) then
!       msg='paw_ij1()%cplex and paw_an1()%cplex must be equal !'
!       MSG_BUG(msg)
!     end if
!     if (pawrhoij_a(1)%cplex<paw_an1(1)%cplex.or.pawrhoij_b(1)%cplex<paw_an1(1)%cplex) then
!       msg='pawrhoij()%cplex must be >=paw_an1()%cplex  !'
!       MSG_BUG(msg)
!     end if
!     if (pawrhoij_a(1)%nspden/=pawrhoij_b(1)%nspden) then
!       msg='pawrhoij_a()%nspden must =pawrhoij_b()%nspden  !'
!       MSG_BUG(msg)
!     end if
!   end if
! end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!!Init contribution to 1st-order (or 2nd-order) energy
! delta_energy(1:2)=zero

!!For some perturbations, nothing else to do
! if (ipert1==natom+1.or.ipert1==natom+10.or.ipert1==natom+11 .or. &
!& ipert2==natom+1.or.ipert2==natom+10.or.ipert2==natom+11) return

!!Various inits
 opt_compch=0;!optvxc=1;optexc=3
 usecore=0;usetcore=0  ! This is true for phonons and Efield pert.
 usexcnhat=maxval(pawtab(1:ntypat)%usexcnhat)
! delta_energy_xc(1:2)=zero;delta_energy_h(1:2)=zero
! dij(1:2)=zero;ro(1:2)=zero

 npts=pawang%angl_size

 d3exc = zero

!================ Loop on atomic sites =======================
 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

   itypat=pawrhoij_1(iatom)%itypat
   mesh_size=pawtab(itypat)%mesh_size
   nspden=pawrhoij_1(iatom)%nspden
   cplex_1=pawrhoij_1(iatom)%cplex
   cplex_2=pawrhoij_2(iatom)%cplex
   cplex_3=pawrhoij_3(iatom)%cplex
   lm_size_all=paw_an0(iatom)%lm_size

!   NOTE : Why lm_size could be dependent of pert?
!   cplex_dijh1=paw_ij1(iatom)%cplex
!   lm_size_1=paw_an1(iatom)%lm_size
!   if (ipert2<=0) lm_size_b=paw_an0(iatom)%lm_size
!   if (ipert2> 0) lm_size_b=paw_an1(iatom)%lm_size

!!  If Vxc potentials are not in memory, compute them
!   if (paw_an1(iatom)%has_vxc/=2) then
   ABI_ALLOCATE(lmselect_tmp,(lm_size_all))
   lmselect_tmp(:)=.true.
!   if (nzlmopt_a==1) lmselect_tmp(:)=lmselect_a(:)

!  Compute on-site 1st-order densities (pert1)
   ABI_ALLOCATE(lmselect_1,(lm_size_all))
   lmselect_1(:)=paw_an0(iatom)%lmselect(:)
   ABI_ALLOCATE(rho1_1,(cplex_1*mesh_size,lm_size_all,nspden))
   ABI_ALLOCATE(trho1_1,(cplex_1*mesh_size,lm_size_all,nspden))
   ABI_ALLOCATE(nhat1_1,(cplex_1*mesh_size,lm_size_all,nspden*usexcnhat))
   call pawdensities(compch,cplex_1,iatom_tot,lmselect_tmp,lmselect_1,&
&   lm_size_all,nhat1_1,nspden,nzlmopt,opt_compch,1-usexcnhat,-1,0,pawang,pawprtvol,&
&   pawrad(itypat),pawrhoij_1(iatom),pawtab(itypat),rho1_1,trho1_1)
!  Compute on-site 1st-order densities (pert2)
   ABI_ALLOCATE(lmselect_2,(lm_size_all))
   lmselect_2(:)=paw_an0(iatom)%lmselect(:)
   ABI_ALLOCATE(rho1_2,(cplex_2*mesh_size,lm_size_all,nspden))
   ABI_ALLOCATE(trho1_2,(cplex_2*mesh_size,lm_size_all,nspden))
   ABI_ALLOCATE(nhat1_2,(cplex_2*mesh_size,lm_size_all,nspden*usexcnhat))
   call pawdensities(compch,cplex_2,iatom_tot,lmselect_tmp,lmselect_2,&
&   lm_size_all,nhat1_2,nspden,nzlmopt,opt_compch,1-usexcnhat,-1,0,pawang,pawprtvol,&
&   pawrad(itypat),pawrhoij_2(iatom),pawtab(itypat),rho1_2,trho1_2)
!  Compute on-site 1st-order densities (pert3)
   ABI_ALLOCATE(lmselect_3,(lm_size_all))
   lmselect_3(:)=paw_an0(iatom)%lmselect(:)
   ABI_ALLOCATE(rho1_3,(cplex_3*mesh_size,lm_size_all,nspden))
   ABI_ALLOCATE(trho1_3,(cplex_3*mesh_size,lm_size_all,nspden))
   ABI_ALLOCATE(nhat1_3,(cplex_3*mesh_size,lm_size_all,nspden*usexcnhat))
   call pawdensities(compch,cplex_3,iatom_tot,lmselect_tmp,lmselect_3,&
&   lm_size_all,nhat1_3,nspden,nzlmopt,opt_compch,1-usexcnhat,-1,0,pawang,pawprtvol,&
&   pawrad(itypat),pawrhoij_3(iatom),pawtab(itypat),rho1_3,trho1_3)

   ABI_DEALLOCATE(lmselect_tmp)

   call paw_dfptnl_xc(cplex_1,cplex_2,cplex_3,d3exc1_iat,ixc,paw_an0(iatom)%k3xc1,lm_size_all,&
&                 lmselect_1,lmselect_2,lmselect_3,nhat1_1,nhat1_2,nhat1_3,&
&                 paw_an0(iatom)%nk3xc1,mesh_size,nspden,pawang,pawrad(itypat),&
                  rho1_1,rho1_2,rho1_3,0)
   d3exc = d3exc + d3exc1_iat

   call paw_dfptnl_xc(cplex_1,cplex_2,cplex_3,d3exc1_iat,ixc,paw_an0(iatom)%k3xc1,lm_size_all,&
&                 lmselect_1,lmselect_2,lmselect_3,nhat1_1,nhat1_2,nhat1_3,&
&                 paw_an0(iatom)%nk3xc1,mesh_size,nspden,pawang,pawrad(itypat),&
                  trho1_1,trho1_2,rho1_3,usexcnhat)
   d3exc = d3exc - d3exc1_iat

!   lm_size_eff=min(lm_size_all,pawang%ylm_size)

!   ABI_ALLOCATE(rho1arr,(cplex_1*nrad,nspden))

!!  Do loop on the angular part (theta,phi)
!   do ipts=1,npts

!  !  Copy the input 1st-order density for this (theta,phi)
!     rho1arr(:,:)=zero
!     if (usexcnhat==0.or.usexcnhat==1) then
!       do ispden=1,nspden
!         do ilm=1,lm_size_eff
!           if (lmselect_1(ilm)) rho1arr(:,ispden)=rho1arr(:,ispden) &
!  &         +rhor1(:,ilm,ispden)*pawang%ylmr(ilm,ipts)
!         end do
!       end do
!     else
!       do ispden=1,nspden
!         do ilm=1,lm_size_eff
!           if (lmselect_1(ilm)) rho1arr(:,ispden)=rho1arr(:,ispden) &
!  &         +(rhor1(:,ilm,ispden)+nhat1(:,ilm,ispden))*pawang%ylmr(ilm,ipts)
!         end do
!       end do
!     end if

!!     if (usecore==1) then
!!       rho1arr(:,1)=rho1arr(:,1)+pawtab(itypat)%coredens(:)
!!       if (nspden==2) rho1arr(:,2)=rho1arr(:,2)+half*corexc1(:)
!!     end if

!  !  ----------------------------------------------------------------------
!  !  ----- Accumulate and store 1st-order change of XC potential
!  !  ----------------------------------------------------------------------

!     if (option/=3) then
!  !    Non-spin-polarized
!       if(nspden==1)then
!         if (cplex_vxc==1) then
!           if (cplex_den==1) then  ! cplex_vxc==1 and cplex_den==1
!             vxc1_(1:nrad,1)=kxc(1:nrad,ipts,1)*rho1arr(1:nrad,1)
!           else                    ! cplex_vxc==1 and cplex_den==2
!             do ir=1,nrad
!               vxc1_(ir,1)=kxc(ir,ipts,1)*rho1arr(2*ir-1,1)
!             end do
!           end if
!         else
!           if (cplex_den==1) then  ! cplex_vxc==2 and cplex_den==1
!             do ir=1,nrad
!               vxc1_(2*ir-1,1)=kxc(ir,ipts,1)*rho1arr(ir,1)
!               vxc1_(2*ir  ,1)=zero
!             end do
!           else                    ! cplex_vxc==2 and cplex_den==2
!             do ir=1,nrad
!               vxc1_(2*ir-1,1)=kxc(ir,ipts,1)*rho1arr(2*ir-1,1)
!               vxc1_(2*ir  ,1)=kxc(ir,ipts,1)*rho1arr(2*ir  ,1)
!             end do
!           end if
!         end if
!  !      Spin-polarized
!       else
!         if (cplex_vxc==1) then
!           if (cplex_den==1) then  ! cplex_vxc==1 and cplex_den==1
!             do ir=1,nrad
!               rho_up=rho1arr(ir,2);rho_dn=rho1arr(ir,1)-rho_up
!               vxc1_(ir,1)=kxc(ir,ipts,1)*rho_up+kxc(ir,ipts,2)*rho_dn
!               vxc1_(ir,2)=kxc(ir,ipts,2)*rho_up+kxc(ir,ipts,3)*rho_dn
!             end do
!           else                    ! cplex_vxc==1 and cplex_den==2
!             do ir=1,nrad
!               jr=2*ir-1
!               rho_up=rho1arr(jr,2);rho_dn=rho1arr(jr,1)-rho_up
!               vxc1_(ir,1)=kxc(ir,ipts,1)*rho_up+kxc(ir,ipts,2)*rho_dn
!               vxc1_(ir,2)=kxc(ir,ipts,2)*rho_up+kxc(ir,ipts,3)*rho_dn
!             end do
!           end if
!         else
!           if (cplex_den==1) then  ! cplex_vxc==2 and cplex_den==1
!             do ir=1,nrad
!               jr=2*ir-1
!               rho_up=rho1arr(ir,2);rho_dn=rho1arr(ir,1)-rho_up
!               vxc1_(jr,1)=kxc(ir,ipts,1)*rho_up+kxc(ir,ipts,2)*rho_dn
!               vxc1_(jr,2)=kxc(ir,ipts,2)*rho_up+kxc(ir,ipts,3)*rho_dn
!             end do
!           else                    ! cplex_vxc==2 and cplex_den==2
!             do ir=1,nrad
!               jr=2*ir
!               rho_up  =rho1arr(jr-1,2);rho_dn  =rho1arr(jr-1,1)-rho_up
!               rhoim_up=rho1arr(jr  ,2);rhoim_dn=rho1arr(jr  ,1)-rhoim_up
!               vxc1_(jr-1,1)=kxc(ir,ipts,1)*rho_up  +kxc(ir,ipts,2)*rho_dn
!               vxc1_(jr  ,1)=kxc(ir,ipts,1)*rhoim_up+kxc(ir,ipts,2)*rhoim_dn
!               vxc1_(jr-1,2)=kxc(ir,ipts,2)*rho_up  +kxc(ir,ipts,3)*rho_dn
!               vxc1_(jr  ,2)=kxc(ir,ipts,2)*rhoim_up+kxc(ir,ipts,3)*rhoim_dn
!             end do
!           end if
!         end if
!       end if

!       if (option<=1) then
!         vxc1(1:cplex_vxc*nrad,ipts,1:nspden)=vxc1_(1:cplex_vxc*nrad,1:nspden)
!       end if

!     else  ! option==3
!       vxc1_(1:cplex_vxc*nrad,1:nspden)=vxc1(1:cplex_vxc*nrad,ipts,1:nspden)
!     end if

!  !  ----------------------------------------------------------------------
!  !  ----- Accumulate and store 2nd-order change of XC energy
!  !  ----------------------------------------------------------------------
!     if (option/=1) then

!  !    For usexnhat=1 particular case, add now compensation density
!       if (usexcnhat==1) then
!         do ispden=1,nspden
!           do ilm=1,lm_size_eff
!             if (lmselect(ilm)) rho1arr(:,ispden)=rho1arr(:,ispden)+nhat1(:,ilm,ispden)*pawang%ylmr(ilm,ipts)
!           end do
!         end do
!       end if

!  !    ----- Calculate d2Exc=Int[Vxc^(1)^*(r).n^(1)(r).dr]
!       LIBPAW_ALLOCATE(ff,(nrad))
!       if (need_impart) then
!         LIBPAW_ALLOCATE(gg,(nrad))
!       end if

!  !    COLLINEAR MAGNETISM
!       if (nspden/=4) then
!         if (cplex_vxc==1.and.cplex_den==1) then       ! cplex_vxc==1 and cplex_den==1
!           ff(:)=vxc1_(:,1)*rho1arr(:,nspden)
!           if (nspden==2) ff(:)=ff(:)+vxc1_(:,2)*(rho1arr(:,1)-rho1arr(:,2))
!           if (need_impart) gg(:)=zero
!         else if (cplex_vxc==2.and.cplex_den==2) then  ! cplex_vxc==2 and cplex_den==2
!           if (.not.need_impart) then      ! Real part only
!             do ir=1,nrad
!               jr=2*ir;v11r=vxc1_(jr-1,1);v11i=vxc1_(jr,1)
!               ro11r=rho1arr(jr-1,nspden);ro11i=rho1arr(jr,nspden)
!               ff(ir)=v11r*ro11r+v11i*ro11i
!             end do
!             if (nspden==2) then
!               do ir=1,nrad
!                 jr=2*ir;v22r=vxc1_(jr-1,2);v22i=vxc1_(jr,2)
!                 ro22r=rho1arr(jr-1,1)-rho1arr(jr-1,2)
!                 ro22i=rho1arr(jr  ,1)-rho1arr(jr  ,2)
!                 ff(ir)=ff(ir)+v22r*ro22r+v22i*ro22i
!               end do
!             end if
!           else
!             do ir=1,nrad                  ! Real and imaginary parts
!               jr=2*ir;v11r=vxc1_(jr-1,1);v11i=vxc1_(jr,1)
!               ro11r=rho1arr(jr-1,nspden);ro11i=rho1arr(jr,nspden)
!               ff(ir)=v11r*ro11r+v11i*ro11i
!               gg(ir)=v11r*ro11i-v11i*ro11r
!             end do
!             if (nspden==2) then
!               do ir=1,nrad
!                 jr=2*ir;v22r=vxc1_(jr-1,2);v22i=vxc1_(jr,2)
!                 ro22r=rho1arr(jr-1,1)-rho1arr(jr-1,2)
!                 ro22i=rho1arr(jr  ,1)-rho1arr(jr  ,2)
!                 ff(ir)=ff(ir)+v22r*ro22r+v22i*ro22i
!                 gg(ir)=gg(ir)+v22r*ro22i-v22i*ro22r
!               end do
!             end if
!           end if
!         else                                          ! other cases for cplex_vxc and cplex_den
!           v11i=zero;ro11i=zero
!           do ir=1,nrad
!             jr=cplex_vxc*(ir-1)+1
!             v11r=vxc1_(jr,1);if (cplex_vxc==2) v11i=vxc1_(jr+1,1)
!             jr=cplex_den*(ir-1)+1
!             ro11r=rho1arr(jr,nspden);if (cplex_den==2) ro11i=rho1arr(jr+1,nspden)
!             ff(ir)=v11r*ro11r+v11i*ro11i
!             if (need_impart) gg(ir)=v11r*ro11i-v11i*ro11r
!           end do
!           if (nspden==2) then
!             v22i=zero;ro22i=zero
!             do ir=1,nrad
!               jr=cplex_vxc*(ir-1)+1
!               v22r=vxc1_(jr,2);if (cplex_vxc==2) v22i=vxc1_(jr+1,2)
!               jr=cplex_den*(ir-1)+1
!               ro22r=rho1arr(jr,1)-rho1arr(jr,2)
!               if (cplex_den==2) ro22i=rho1arr(jr+1,1)-rho1arr(jr+1,2)
!               ff(ir)=ff(ir)+v22r*ro22r+v22i*ro22i
!               gg(ir)=gg(ir)+v22r*ro22i-v22i*ro22r
!             end do
!           end if
!         end if ! cplex_vxc and cplex_den

!!  !      NON-COLLINEAR MAGNETISM
!!       else
!!         if (cplex_vxc==1.and.cplex_den==1) then   ! cplex_vxc==1 and cplex_den==1
!!           ff(:)=half*(vxc1_(:,1)*(rho1arr(:,1)+rho1arr(:,4)) &
!!  &         +vxc1_(:,2)*(rho1arr(:,1)-rho1arr(:,4))) &
!!  &         +vxc1_(:,3)*rho1arr(:,2) &
!!  &         -vxc1_(:,4)*rho1arr(:,3)
!!           if (need_impart) gg(:)=zero
!!         else                                      ! other cases for cplex_vxc and cplex_den

!!  !        V is stored as : v^11, v^22, V^12, i.V^21 (each are complex)
!!  !        N is stored as : n, m_x, m_y, mZ          (each are complex)
!!           do ir=1,nrad
!!             jr=cplex_vxc*(ir-1)+1
!!             v11r= vxc1_(jr,1);v22r= vxc1_(jr,2)
!!             v12r= vxc1_(jr,3);v21i=-vxc1_(jr,1)
!!             if (cplex_vxc==2) then
!!               v11i= vxc1_(jr+1,1);v22i= vxc1_(jr+1,2)
!!               v12i= vxc1_(jr+1,3);v21r= vxc1_(jr+1,1)
!!             else
!!               v11i=zero;v22i=zero
!!               v12i=zero;v21i=zero
!!             end if
!!             jr=cplex_den*(ir-1)+1
!!             ro11r= rho1arr(jr,1)+rho1arr(jr,4)
!!             ro22r= rho1arr(jr,1)-rho1arr(jr,4)
!!             ro12r= rho1arr(jr,2);ro12i=-rho1arr(jr,3)
!!             ro21r= rho1arr(jr,2);ro21i= rho1arr(jr,3)
!!             if (cplex_den==2) then
!!               ro11i=rho1arr(jr+1,1)+rho1arr(jr+1,4)
!!               ro22i=rho1arr(jr+1,1)-rho1arr(jr+1,4)
!!               ro12r=ro12r+rho1arr(jr+1,3);ro12i=ro12i+rho1arr(jr+1,2)
!!               ro21r=ro21r-rho1arr(jr+1,3);ro21i=ro21i+rho1arr(jr+1,2)
!!             else
!!               ro11i=zero;ro22i=zero
!!             end if
!!  !          Real part
!!             ff(ir)=half*(v11r*ro11r+v11i*ro11i+v22r*ro22r+v22i*ro22i &
!!  &           +v12r*ro12r+v12i*ro12i+v21r*ro21r+v21i*ro21i)
!!  !          Imaginary part
!!             if (need_impart) gg(ir)=half*(v11r*ro11i-v11i*ro11r+v22r*ro22i-v22i*ro22r &
!!  &           +v12r*ro12i-v12i*ro12r+v21r*ro21i-v21i*ro21r)
!!           end do
!!         end if ! cplex_vxc and cplex_den
!!       end if ! nspden

!       ff(:)=ff(:)*pawrad%rad(:)**2
!       call simp_gen(vxcrho,ff,pawrad)
!       d2enxc=d2enxc+vxcrho*pawang%angwgth(ipts)
!       LIBPAW_DEALLOCATE(ff)

!       if (need_impart) then
!         gg(:)=gg(:)*pawrad%rad(:)**2
!         call simp_gen(vxcrho,gg,pawrad)
!         d2enxc_im=d2enxc_im+vxcrho*pawang%angwgth(ipts)
!         LIBPAW_DEALLOCATE(gg)
!       end if

!     end if

!!  ----- End of the loop on npts (angular part)
!   end do

!!  Add the four*pi factor of the angular integration
!   d3exc=d3exc*four_pi

   ABI_DEALLOCATE(lmselect_1)
   ABI_DEALLOCATE(lmselect_2)
   ABI_DEALLOCATE(lmselect_3)
   ABI_DEALLOCATE(nhat1_1)
   ABI_DEALLOCATE(nhat1_2)
   ABI_DEALLOCATE(nhat1_3)
   ABI_DEALLOCATE(rho1_1)
   ABI_DEALLOCATE(rho1_2)
   ABI_DEALLOCATE(rho1_3)
   ABI_DEALLOCATE(trho1_1)
   ABI_DEALLOCATE(trho1_2)
   ABI_DEALLOCATE(trho1_3)

!!  Compute contribution to 1st-order(or 2nd-order) energy from 1st-order Hartree potential
!   nspdiag=1;if (nspden==2) nspdiag=2
!   do ispden=1,nspdiag
!     if (cplex_dijh1==1) then
!       jrhoij=1
!       do irhoij=1,pawrhoij_b(iatom)%nrhoijsel
!         klmn=pawrhoij_b(iatom)%rhoijselect(irhoij)
!         dij(1)=paw_ij1(iatom)%dijhartree(klmn)
!         ro(1)=pawrhoij_b(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
!         delta_energy_h(1)=delta_energy_h(1)+ro(1)*dij(1)
!         if (cplex_b==2) then
!           ro(2)=pawrhoij_b(iatom)%rhoijp(jrhoij+1,ispden)*pawtab(itypat)%dltij(klmn)
!           delta_energy_h(2)=delta_energy_h(2)+ro(2)*dij(1)
!         end if
!         jrhoij=jrhoij+cplex_b
!       end do
!     else ! cplex_dijh1==2
!       jrhoij=1
!       do irhoij=1,pawrhoij_b(iatom)%nrhoijsel
!         klmn=pawrhoij_b(iatom)%rhoijselect(irhoij)
!         dij(1:2)=paw_ij1(iatom)%dijhartree(2*klmn-1:2*klmn)
!         ro(1)=pawrhoij_b(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
!         delta_energy_h(1)=delta_energy_h(1)+ro(1)*dij(1)
!         delta_energy_h(2)=delta_energy_h(2)-ro(1)*dij(2)
!         if (cplex_b==2) then
!           ro(2)=pawrhoij_b(iatom)%rhoijp(jrhoij+1,ispden)*pawtab(itypat)%dltij(klmn)
!           delta_energy_h(1)=delta_energy_h(1)+ro(2)*dij(2)
!           delta_energy_h(2)=delta_energy_h(2)+ro(2)*dij(1)
!         end if
!         jrhoij=jrhoij+cplex_b
!       end do
!     end if
!   end do

!  ================ End loop oon atomic sites =======================
 end do

!!Final building of 1st-order (or 2nd-order) energy
! delta_energy(1:2)=delta_energy_xc(1:2)+delta_energy_h(1:2)

!!Reduction in case of parallelism
! if (paral_atom) then
!   call xmpi_sum(delta_energy,my_comm_atom,ierr)
! end if

!!Destroy atom table used for parallelism
! call free_my_atmtab(my_atmtab,my_atmtab_allocated)

! call timab(567,2,tsec)

 DBG_EXIT("COLL")

end subroutine paw_dfptnl_energy
!!***

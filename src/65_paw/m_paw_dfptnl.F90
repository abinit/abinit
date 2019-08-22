!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_dfptnl
!! NAME
!!  m_paw_dfptnl
!!
!! FUNCTION
!!  This module contains several routines used to compute PAW contributions to a 3rd-order energy
!!   or 2nd-order PAW occupancies.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2019 ABINIT group (LB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_dfptnl

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi

 use m_pawang,     only : pawang_type
 use m_pawrad,     only : pawrad_type,simp_gen
 use m_pawtab,     only : pawtab_type
 use m_paw_an,     only : paw_an_type
 use m_pawrhoij,   only : pawrhoij_type
 use m_pawcprj,    only : pawcprj_type
 use m_paw_denpot, only : pawdensities
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab

 implicit none

 private

!public procedures.
 public :: paw_dfptnl_energy   ! Compute the XC PAW on-site contributions to a 3rd-order energy
 public :: paw_dfptnl_xc       ! Compute a contribution of the 3rd-derivative of XC energy of ONE PAW sphere
 public :: paw_dfptnl_accrhoij ! Accumulate the 2nd order PAW quantities rhoij^(2)

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_dfptnl/paw_dfptnl_energy
!! NAME
!! paw_dfptnl_energy
!!
!! FUNCTION
!! Compute the XC PAW on-site contributions to a 3rd-order energy.
!! It is equal to:
!!    E_onsite= \sum_at [ E_at(kxc,rho1,rho2,rho3) - E_at(tkxc,trho1,trho2,trho3) ]
!! where E_at(...) is computed in paw_dfptnl_xc.F90.
!! The atomic densities are computed from pawrhoij1,pawrhoij2 and pawrhoij3.
!! This routine is similar to pawdfptenergy.F90 but is implemented independently
!! in order to not overload the original routine.
!! LDA ONLY - USE THE DENSITY OVER A WHOLE SPHERICAL GRID (r,theta,phi)
!!
!! INPUTS
!!  ixc= choice of exchange-correlation scheme
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  ntypat=number of types of atoms in unit cell.
!!  paw_an0(natom) <type(paw_an_type)>=paw arrays for 0th-order quantities given on angular mesh
!!  paw_an1(natom) <type(paw_an_type)>=paw arrays for 1st-order quantities given on angular mesh
!!                                     This corresponds to (j1) perturbation
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij_1-2-3(natom) <type(pawrhoij_type)>= paw rhoij 1st-order occupancies
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!
!! OUTPUT
!!  d3exc= real and imaginary parts of the contribution to the third derivative of the total energy
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      dfptnl_pert
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,pawdensities,pawxc_dfpt,pawxcm3
!!      timab,xmpi_sum
!!
!! SOURCE

subroutine paw_dfptnl_energy(d3exc,ixc,my_natom,natom,ntypat,&
&                    paw_an0,pawang,pawprtvol,pawrad,&
&                    pawrhoij_1,pawrhoij_2,pawrhoij_3,&
&                    pawtab,pawxcdev,&
&                    mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ixc,my_natom,natom,ntypat
 integer,intent(in) :: pawprtvol,pawxcdev
 integer,optional,intent(in) :: comm_atom
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(out) :: d3exc(2)
 type(paw_an_type),intent(in) :: paw_an0(my_natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij_1(my_natom)
 type(pawrhoij_type),intent(in) :: pawrhoij_2(my_natom)
 type(pawrhoij_type),intent(in) :: pawrhoij_3(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: cplex_1,cplex_2,cplex_3,iatom,iatom_tot,itypat
 integer :: lm_size_all,mesh_size,my_comm_atom,npts,nspden,nzlmopt
 integer :: opt_compch,usecore,usetcore,usexcnhat
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: compch,d3exc1_iat(2)
 character(len=500) :: msg
!arrays
 integer,pointer :: my_atmtab(:)
 logical,allocatable :: lmselect_1(:),lmselect_2(:),lmselect_3(:),lmselect_tmp(:)
 real(dp),allocatable :: nhat1_1(:,:,:),rho1_1(:,:,:),trho1_1(:,:,:)
 real(dp),allocatable :: nhat1_2(:,:,:),rho1_2(:,:,:),trho1_2(:,:,:)
 real(dp),allocatable :: nhat1_3(:,:,:),rho1_3(:,:,:),trho1_3(:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 nzlmopt = 0 ! compute all LM-moments of the density and use all LM-moments

 if (pawxcdev/=0) then
   msg="paw_dfptnl_energy is not implemented for pawxcdev/=0"
   MSG_BUG(msg)
 end if
 if (my_natom>0) then
   if (pawrhoij_1(1)%qphase/=1.or.pawrhoij_2(1)%qphase/=1.or.pawrhoij_3(1)%qphase/=1) then
     msg="paw_dfptnl_energy not supposed to be called with q/=0!"
     MSG_BUG(msg)
   end if
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!!Various inits
 opt_compch=0; !optvxc=1;optexc=3
 usecore=0;usetcore=0  ! This is true for phonons and Efield pert.
 usexcnhat=maxval(pawtab(1:ntypat)%usexcnhat)

 npts=pawang%angl_size

 d3exc = zero

!================ Loop on atomic sites =======================
 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

   itypat=pawrhoij_1(iatom)%itypat
   mesh_size=pawtab(itypat)%mesh_size
   nspden=pawrhoij_1(iatom)%nspden
   cplex_1=pawrhoij_1(iatom)%cplex_rhoij
   cplex_2=pawrhoij_2(iatom)%cplex_rhoij
   cplex_3=pawrhoij_3(iatom)%cplex_rhoij
   lm_size_all=paw_an0(iatom)%lm_size

   ABI_ALLOCATE(lmselect_tmp,(lm_size_all))
   lmselect_tmp(:)=.true.

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
&   lmselect_1,lmselect_2,lmselect_3,nhat1_1,nhat1_2,nhat1_3,&
&   paw_an0(iatom)%nk3xc1,mesh_size,nspden,pawang,pawrad(itypat),&
&   rho1_1,rho1_2,rho1_3,0)
   d3exc = d3exc + d3exc1_iat

   call paw_dfptnl_xc(cplex_1,cplex_2,cplex_3,d3exc1_iat,ixc,paw_an0(iatom)%k3xct1,lm_size_all,&
&   lmselect_1,lmselect_2,lmselect_3,nhat1_1,nhat1_2,nhat1_3,&
&   paw_an0(iatom)%nk3xc1,mesh_size,nspden,pawang,pawrad(itypat),&
&   trho1_1,trho1_2,trho1_3,usexcnhat)
   d3exc = d3exc - d3exc1_iat

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

!  ================ End loop oon atomic sites =======================
 end do

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

!----------------------------------------------------------------------

!!****f* m_paw_dfptnl/paw_dfptnl_xc
!! NAME
!! paw_dfptnl_xc
!!
!! FUNCTION
!! Compute a contribution of the third derivative of XC energy of ONE PAW sphere Om_a.
!! It is equal to:
!!   E_at(kxc,rho1,rho2,rho3) = Int_{Om_a} dr kxc(r) * rho1(r) * rho2(r) * rho3(r)
!! where kxc,rho1,rho2 and rho3 are inputs.
!! This routine is similar to m_pawxc.F90:pawxc_dfpt(...) but is implemented independently
!! in order to not overload the original routine.
!! LDA ONLY - USE THE DENSITY OVER A WHOLE SPHERICAL GRID (r,theta,phi)
!!
!! INPUTS
!!  cplex_1-2-3= if 1, 1st-order densities are REAL, if 2, COMPLEX
!!  d3exc1_iat=third-order derivative to compute
!!  ixc= choice of exchange-correlation scheme
!!  kxc(nrad,pawang%angl_size,nkxc)=GS xc kernel
!!  lm_size=size of density array rhor (see below)
!!  lmselect1-2-3(lm_size)=select the non-zero LM-moments of input density rhor1-2-3
!!  nhat1-2-3(cplex_den*nrad,lm_size,nspden)=first-order change of compensation density
!!                                        (total in 1st half and spin-up in 2nd half if nspden=2)
!!  nkxc=second dimension of the kxc array
!!  nrad=size of radial mesh for densities/potentials (might be different from pawrad%mesh_size)
!!  nspden=number of spin-density components
!!  option=0  compute both 2nd-order XC energy and 1st-order potential
!!         1  compute only 1st-order XC potential
!!         2  compute only 2nd-order XC energy, XC potential is temporary computed here
!!         3  compute only 2nd-order XC energy, XC potential is input in vxc1(:)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  rhor1-2-3(cplex_den*nrad,lm_size,nspden)=first-order change of density
!!  usexcnhat= 0 if compensation density does not have to be used
!!             1 if compensation density has to be used in d2Exc only
!!
!! OUTPUT
!!  d3exc1_iat = E_at(kxc,rho1,rho2,rho3) (see FUNCTION above)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      paw_dfptnl_energy
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_dfptnl_xc(cplex_1,cplex_2,cplex_3,d3exc1_iat,ixc,kxc,lm_size,lmselect1,lmselect2,lmselect3,&
&                 nhat1,nhat2,nhat3,nkxc,nrad,nspden,pawang,pawrad,rhor1,rhor2,rhor3,usexcnhat)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex_1,cplex_2,cplex_3,ixc,lm_size,nkxc,nrad,nspden,usexcnhat
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect1(lm_size),lmselect2(lm_size),lmselect3(lm_size)
 real(dp),intent(out) :: d3exc1_iat(2)
 real(dp),intent(in) :: kxc(nrad,pawang%angl_size,nkxc)
 real(dp),intent(in) :: nhat1(cplex_1*nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: nhat2(cplex_2*nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: nhat3(cplex_3*nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: rhor1(cplex_1*nrad,lm_size,nspden)
 real(dp),intent(in) :: rhor2(cplex_2*nrad,lm_size,nspden)
 real(dp),intent(in) :: rhor3(cplex_3*nrad,lm_size,nspden)

!Local variables-------------------------------
!scalars
 integer :: ii,ilm,ipts,ispden,lm_size_eff,npts
 real(dp) :: d3exc1_int,rho1u,rho1d,rho2u,rho2d,rho3u,rho3d
 character(len=500) :: msg
!arrays
! real(dp) :: tsec(2)
 real(dp),allocatable :: ff(:),rho1arr(:,:),rho2arr(:,:),rho3arr(:,:)

! *************************************************************************

!----------------------------------------------------------------------
!----- Initializations
!----------------------------------------------------------------------

 npts=pawang%angl_size
 lm_size_eff=min(lm_size,pawang%ylm_size)

 d3exc1_iat(:) = zero

!Special case: no XC applied
 if (ixc==0) then
   msg='Note that no xc is applied (ixc=0). Returning'
   MSG_WARNING(msg)
   return
 end if

 ABI_ALLOCATE(rho1arr,(cplex_1*nrad,nspden))
 ABI_ALLOCATE(rho2arr,(cplex_2*nrad,nspden))
 ABI_ALLOCATE(rho3arr,(cplex_3*nrad,nspden))

!Restriction : all cplex must be 1
 if (cplex_1/=1.or.cplex_2/=1.or.cplex_3/=1) then
   msg='All cplex must be one (for the moment...)'
   MSG_BUG(msg)
 end if
!Restriction : nspden must be 1
! if (nkxc>1) then
!   msg='nkxc must be one (<=> nspden=1) (for the moment...)'
!   MSG_BUG(msg)
! end if

 ABI_ALLOCATE(ff,(nrad))

!!----------------------------------------------------------------------
!!----- Loop on the angular part and inits
!!----------------------------------------------------------------------

!Do loop on the angular part (theta,phi)
 do ipts=1,npts

!  Copy the input 1st-order density for this (theta,phi) - PERT1
   rho1arr(:,:)=zero
   if (usexcnhat==0) then
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect1(ilm)) rho1arr(:,ispden)=rho1arr(:,ispden) &
&         +rhor1(:,ilm,ispden)*pawang%ylmr(ilm,ipts)
       end do
     end do
   else
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect1(ilm)) rho1arr(:,ispden)=rho1arr(:,ispden) &
&         +(rhor1(:,ilm,ispden)+nhat1(:,ilm,ispden))*pawang%ylmr(ilm,ipts)
       end do
     end do
   end if

!  Copy the input 1st-order density for this (theta,phi) - PERT2
   rho2arr(:,:)=zero
   if (usexcnhat==0) then
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect2(ilm)) rho2arr(:,ispden)=rho2arr(:,ispden) &
&         +rhor2(:,ilm,ispden)*pawang%ylmr(ilm,ipts)
       end do
     end do
   else
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect2(ilm)) rho2arr(:,ispden)=rho2arr(:,ispden) &
&         +(rhor2(:,ilm,ispden)+nhat2(:,ilm,ispden))*pawang%ylmr(ilm,ipts)
       end do
     end do
   end if

!  Copy the input 1st-order density for this (theta,phi) - PERT3
   rho3arr(:,:)=zero
   if (usexcnhat==0) then
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect3(ilm)) rho3arr(:,ispden)=rho3arr(:,ispden) &
&         +rhor3(:,ilm,ispden)*pawang%ylmr(ilm,ipts)
       end do
     end do
   else
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect3(ilm)) rho3arr(:,ispden)=rho3arr(:,ispden) &
&         +(rhor3(:,ilm,ispden)+nhat3(:,ilm,ispden))*pawang%ylmr(ilm,ipts)
       end do
     end do
   end if

!  ----------------------------------------------------------------------
!  ----- Accumulate and store 3nd-order change of XC energy
!  ----------------------------------------------------------------------

   if (cplex_1==1.and.cplex_2==1.and.cplex_3==1) then ! all cplex are 1 :
     if (nspden==1) then
       ff(:)=kxc(:,ipts,1)*rho1arr(:,1)*rho2arr(:,1)*rho3arr(:,1)
     else if (nspden==2) then
       do ii=1,nrad
         rho1u=rho1arr(ii,2)
         rho1d=rho1arr(ii,1)-rho1arr(ii,2)
         rho2u=rho2arr(ii,2)
         rho2d=rho2arr(ii,1)-rho2arr(ii,2)
         rho3u=rho3arr(ii,2)
         rho3d=rho3arr(ii,1)-rho3arr(ii,2)
         ff(ii)=&
!          uuu                                uud
&         kxc(ii,ipts,1)*rho1u*rho2u*rho3u + kxc(ii,ipts,2)*rho1u*rho2u*rho3d + &
!          udu                                udd
&         kxc(ii,ipts,2)*rho1u*rho2d*rho3u + kxc(ii,ipts,3)*rho1u*rho2d*rho3d + &
!          duu                                dud
&         kxc(ii,ipts,2)*rho1d*rho2u*rho3u + kxc(ii,ipts,3)*rho1d*rho2u*rho3d + &
!          ddu                                ddd
&         kxc(ii,ipts,3)*rho1d*rho2d*rho3u + kxc(ii,ipts,4)*rho1d*rho2d*rho3d
       end do
     end if
   end if

   ff(1:nrad)=ff(1:nrad)*pawrad%rad(1:nrad)**2
   call simp_gen(d3exc1_int,ff,pawrad)
   d3exc1_iat(1)=d3exc1_iat(1)+d3exc1_int*pawang%angwgth(ipts)

!  ----- End of the loop on npts (angular part)
 end do

 d3exc1_iat = d3exc1_iat*four_pi

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(rho1arr)
 ABI_DEALLOCATE(rho2arr)
 ABI_DEALLOCATE(rho3arr)

end subroutine paw_dfptnl_xc
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/paw_dfptnl_accrhoij
!!
!! NAME
!! paw_dfptnl_accrhoij
!!
!! FUNCTION
!! Accumulate the 2nd order PAW quantities rhoij^(2) (augmentation occupancies)
!! This routine is similar to pawaccrhoij.F90 but is implemented independently
!! in order to not overload the original routine.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (sorted-->random), inverse of atindx.
!!  cplex: if 1, WFs (or 1st-order WFs) are REAL, if 2, COMPLEX
!!  cwaveprj0_pert1(natom,nspinor) = wave function at given n,k projected with non-local projectors:
!!                                  cwaveprj0%cp    =<p_i|
!!                                  cwaveprj0%dcp(1)=<p_i^(pert1)|Cnk>
!!  cwaveprj0_pert2(natom,nspinor) = wave function at given n,k projected with non-local projectors:
!!                                  cwaveprj0%cp    =<p_i|Cnk>
!!                                  cwaveprj0%dcp(1)=<p_i^(pert1)|Cnk>
!!  cwaveprj1_pert12(natom,nspinor)= 1st order wave function at given n,k projected with non-local projectors:
!!                                  cwaveprj1%cp    =<p_i|Cnk^(pert1)>
!!                                  cwaveprj1%dcp(1)=<p_i^(pert2)|Cnk^(pert1)>
!!  cwaveprj1_pert21(natom,nspinor)= 1st order wave function at given n,k projected with non-local projectors:
!!                                  cwaveprj1%cp    =<p_i|Cnk^(pert2)>
!!                                  cwaveprj1%dcp(1)=<p_i^(pert1)|Cnk^(pert2)>
!!  ipert1=index of the first perturbation
!!  ipert2=index of the second perturbation
!!  isppol=index of current spin component
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell
!!  nspinor=number of spinorial components (on current proc)
!!  occ_k=occupation number for current band n,k
!!  wtk_k=weight assigned to current k-point
!!
!! SIDE EFFECTS
!!  pawrhoij(natom) <type(pawrhoij_type)>= 2-nd order paw rhoij occupancies and related data
!!  On output, has been updated with the contribution of current n,k
!!        pawrhoij(:)%rhoij_(lmn2_size,nspden) (non symetrized)
!!
!! PARENTS
!!      paw_dfptnl_pert
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab
!!
!! SOURCE

 subroutine paw_dfptnl_accrhoij(atindx,cplex,cwaveprj0_pert1,cwaveprj0_pert2,&
&                       cwaveprj1_pert12,cwaveprj1_pert21,ipert1,ipert2,isppol,my_natom,natom,&
&                       nspinor,occ_k,pawrhoij,wtk_k,&
&                       comm_atom,mpi_atmtab ) ! optional (parallelism)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,ipert1,ipert2,isppol,my_natom,natom,nspinor
 integer,optional,intent(in) :: comm_atom
 real(dp),intent(in) :: occ_k,wtk_k
!arrays
 integer,intent(in) :: atindx(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(pawcprj_type),intent(in) :: cwaveprj0_pert1(natom,nspinor),cwaveprj0_pert2(natom,nspinor)
 type(pawcprj_type),intent(in) :: cwaveprj1_pert12(natom,nspinor),cwaveprj1_pert21(natom,nspinor)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)

!Local variables ---------------------------------------
!scalars
 integer :: cplex_rhoij,iatm,iatom,iatom1,ilmn,iplex,j0lmn,jlmn,klmn,klmn_im,klmn_re
 integer :: my_comm_atom,ncpgr
 logical :: compute_impart,compute_impart_cplex
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: ro11_im,ro11_re,weight
 character(len=500) :: message
!arrays
 integer,pointer :: my_atmtab(:)
 real(dp) :: cpi0(2,nspinor),d1cpi0(2,nspinor),d2cpi0(2,nspinor)
 real(dp) :: cpj0(2,nspinor),d1cpj0(2,nspinor),d2cpj0(2,nspinor)
 real(dp) :: cpi1(2,nspinor),d2cpi1(2,nspinor)
 real(dp) :: cpj1(2,nspinor),d2cpj1(2,nspinor)
 real(dp) :: cpi2(2,nspinor),d1cpi2(2,nspinor)
 real(dp) :: cpj2(2,nspinor),d1cpj2(2,nspinor)

! ***********************************************************************

 DBG_ENTER("COLL")

 if (my_natom==0) return

 ncpgr=1
 if (ipert1<0.or.ipert1>natom+2.or.ipert2<0.or.ipert2>natom+2) then
   message = 'paw_dfptnl_accrhoij: Necessary conditions on ipert1 or ipert2: 0<=ipert<=natom+2'
   MSG_BUG(message)
 end if
 if (pawrhoij(1)%qphase/=1) then
   message="paw_dfptnl_accrhoij not supposed to be called with q/=0!"
   MSG_BUG(message)
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,&
& my_natom_ref=my_natom)

 weight=wtk_k*occ_k
 if (pawrhoij(1)%nspden==2.and.pawrhoij(1)%nsppol==1.and.nspinor==1) weight=half*weight

!!  ==================================================================
!!  === Accumulate (n,k) contribution to partial 2nd-order rhoij   ===
!!  ==================================================================

 compute_impart=(pawrhoij(1)%cplex_rhoij==2)
 compute_impart_cplex=((pawrhoij(1)%cplex_rhoij==2).and.(cplex==2))

! NOT USED FOR PAWRHO21! => only for PAWRHO2 (full second derivative)
!!Accumulate :   < Psi^(pert1) | p_i^(0) > < p_j^(0) | Psi^(pert2) >
!!             + < Psi^(pert2) | p_i^(0) > < p_j^(0) | Psi^(pert1) >
! if (nspinor==1) then
!   do iatom=1,my_natom
!     iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
!     iatm=atindx(iatom1)
!     cplex_rhoij=pawrhoij(iatom)%cplex_rhoij
!     do jlmn=1,pawrhoij(iatom)%lmn_size
!       j0lmn=jlmn*(jlmn-1)/2
!       cpj1(1:2,1)=cwaveprj1_pert12(iatm,1)%cp(1:2,jlmn)   ! < p_j^(0) | Psi^(pert1) >
!       cpj2(1:2,1)=cwaveprj1_pert21(iatm,1)%cp(1:2,jlmn)   ! < p_j^(0) | Psi^(pert2) >
!       do ilmn=1,jlmn
!         klmn=j0lmn+ilmn
!         klmn_re=cplex_rhoij*(klmn-1)+1
!         cpi1(1:2,1)=cwaveprj1_pert12(iatm,1)%cp(1:2,ilmn) ! < p_i^(0) | Psi^(pert1) >
!         cpi2(1:2,1)=cwaveprj1_pert21(iatm,1)%cp(1:2,ilmn) ! < p_i^(0) | Psi^(pert2) >
!         ro11_re=zero
!         do iplex=1,cplex
!           ro11_re=ro11_re+cpi1(iplex,1)*cpj2(iplex,1)+cpj1(iplex,1)*cpi2(iplex,1)
!         end do
!         pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol)+weight*ro11_re
!         if (compute_impart_cplex) then
!           klmn_im=klmn_re+1
!           ro11_im=        cpi1(1,1)*cpj2(2,1)-cpi1(2,1)*cpj2(1,1)
!           ro11_im=ro11_im+cpj1(1,1)*cpi2(2,1)-cpj1(2,1)*cpi2(1,1)
!           pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol)+weight*ro11_im
!         end if
!       end do
!     end do
!   end do
! else ! nspinor=2
!   MSG_BUG("paw_dfptnl_accrhoij is not implemented for nspinor=2")
! end if

!Accumulate :   < Psi^(pert1) | p_i^(pert2) > < p_j^(0)     | Psi^(0)     >
!             + < Psi^(pert1) | p_i^(0)     > < p_j^(pert2) | Psi^(0)     >
!             + < Psi^(0)     | p_i^(pert2) > < p_j^(0)     | Psi^(pert1) >
!             + < Psi^(0)     | p_i^(0)     > < p_j^(pert2) | Psi^(pert1) >
 if (ipert2>0.and.ipert2<=natom) then
   if (nspinor==1) then
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       if (iatom/=ipert2) cycle ! To move atom "ipert2" does not change projectors of other atoms
       cplex_rhoij=pawrhoij(iatom)%cplex_rhoij
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         cpj0(1:2,1)  =cwaveprj0_pert2 (iatm,1)% cp(1:2,  jlmn)   ! < p_j^(0)     | Psi^(0)     >
         d2cpj0(1:2,1)=cwaveprj0_pert2 (iatm,1)%dcp(1:2,1,jlmn)   ! < p_j^(pert2) | Psi^(0)     >
         cpj1(1:2,1)  =cwaveprj1_pert12(iatm,1)% cp(1:2,  jlmn)   ! < p_j^(0)     | Psi^(pert1) >
         d2cpj1(1:2,1)=cwaveprj1_pert12(iatm,1)%dcp(1:2,1,jlmn)   ! < p_j^(pert2) | Psi^(pert1) >
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           cpi0(1:2,1)  =cwaveprj0_pert2 (iatm,1)% cp(1:2,  ilmn) ! < p_i^(0)     | Psi^(0)     >
           d2cpi0(1:2,1)=cwaveprj0_pert2 (iatm,1)%dcp(1:2,1,ilmn) ! < p_i^(pert2) | Psi^(0)     >
           cpi1(1:2,1)  =cwaveprj1_pert12(iatm,1)% cp(1:2,  ilmn) ! < p_i^(0)     | Psi^(pert1) >
           d2cpi1(1:2,1)=cwaveprj1_pert12(iatm,1)%dcp(1:2,1,ilmn) ! < p_i^(pert2) | Psi^(pert1) >
           ro11_re=zero
           do iplex=1,cplex
             ro11_re=ro11_re+d2cpi1(iplex,1)*  cpj0(iplex,1)
             ro11_re=ro11_re+  cpi1(iplex,1)*d2cpj0(iplex,1)
             ro11_re=ro11_re+d2cpi0(iplex,1)*  cpj1(iplex,1)
             ro11_re=ro11_re+  cpi0(iplex,1)*d2cpj1(iplex,1)
           end do
           pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol)+weight*ro11_re
           if (compute_impart_cplex) then
             klmn_im=klmn_re+1
             ro11_im=        d2cpi1(1,1)*  cpj0(2,1)-d2cpi1(2,1)*  cpj0(1,1)
             ro11_im=ro11_im+  cpi1(1,1)*d2cpj0(2,1)-  cpi1(2,1)*d2cpj0(1,1)
             ro11_im=ro11_im+d2cpi0(1,1)*  cpj1(2,1)-d2cpi0(2,1)*  cpj1(1,1)
             ro11_im=ro11_im+  cpi0(1,1)*d2cpj1(2,1)-  cpi0(2,1)*d2cpj1(1,1)
             pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol)+weight*ro11_im
           end if
         end do
       end do
     end do
   else ! nspinor=2
     MSG_BUG("paw_dfptnl_accrhoij is not implemented for nspinor=2")
   end if
 end if

!Accumulate : < Psi^(pert2) | p_i^(pert1) > < p_j^(0)     | Psi^(0)     >
!           + < Psi^(pert2) | p_i^(0)     > < p_j^(pert1) | Psi^(0)     >
!           + < Psi^(0)     | p_i^(pert1) > < p_j^(0)     | Psi^(pert2) >
!           + < Psi^(0)     | p_i^(0) >     < p_j^(pert1) | Psi^(pert2) >
 if (ipert1>0.and.ipert1<=natom) then
   if (nspinor==1) then
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       cplex_rhoij=pawrhoij(iatom)%cplex_rhoij
       if (iatom/=ipert1) cycle ! To move atom "ipert1" does not change projectors of other atoms
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         cpj0(1:2,1)  =cwaveprj0_pert1 (iatm,1)% cp(1:2,  jlmn)   ! < p_j^(0)     | Psi^(0)     >
         d1cpj0(1:2,1)=cwaveprj0_pert1 (iatm,1)%dcp(1:2,1,jlmn)   ! < p_j^(pert1) | Psi^(0)     >
         cpj2(1:2,1)  =cwaveprj1_pert21(iatm,1)% cp(1:2,  jlmn)   ! < p_j^(0)     | Psi^(pert2) >
         d1cpj2(1:2,1)=cwaveprj1_pert21(iatm,1)%dcp(1:2,1,jlmn)   ! < p_j^(pert1) | Psi^(pert2) >
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           cpi0(1:2,1)  =cwaveprj0_pert1 (iatm,1)% cp(1:2,  ilmn) ! < p_i^(0)     | Psi^(0)     >
           d1cpi0(1:2,1)=cwaveprj0_pert1 (iatm,1)%dcp(1:2,1,ilmn) ! < p_i^(pert1) | Psi^(0)     >
           cpi2(1:2,1)  =cwaveprj1_pert21(iatm,1)% cp(1:2,  ilmn) ! < p_i^(0)     | Psi^(pert2) >
           d1cpi2(1:2,1)=cwaveprj1_pert21(iatm,1)%dcp(1:2,1,ilmn) ! < p_i^(pert1) | Psi^(pert2) >
           ro11_re=zero
           do iplex=1,cplex
             ro11_re=ro11_re+d1cpi2(iplex,1)*  cpj0(iplex,1)
             ro11_re=ro11_re+  cpi2(iplex,1)*d1cpj0(iplex,1)
             ro11_re=ro11_re+d1cpi0(iplex,1)*  cpj2(iplex,1)
             ro11_re=ro11_re+  cpi0(iplex,1)*d1cpj2(iplex,1)
           end do
           pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol)+weight*ro11_re
           if (compute_impart_cplex) then
             klmn_im=klmn_re+1
             ro11_im=        d1cpi2(1,1)*  cpj0(2,1)-d1cpi2(2,1)*  cpj0(1,1)
             ro11_im=ro11_im+  cpi2(1,1)*d1cpj0(2,1)-  cpi2(2,1)*d1cpj0(1,1)
             ro11_im=ro11_im+d1cpi0(1,1)*  cpj2(2,1)-d1cpi0(2,1)*  cpj2(1,1)
             ro11_im=ro11_im+  cpi0(1,1)*d1cpj2(2,1)-  cpi0(2,1)*d1cpj2(1,1)
             pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol)+weight*ro11_im
           end if
         end do
       end do
     end do
   else ! nspinor=2
     message="paw_dfptnl_accrhoij is not implemented for nspinor=2!"
     MSG_BUG(message)
   end if
 end if
!  End

!Accumulate :   < Psi^(0) | p_i^(pert1) > < p_j^(pert2) | Psi^(0) >
!             + < Psi^(0) | p_i^(pert2) > < p_j^(pert1) | Psi^(0) >
 if (ipert1>0.and.ipert1<=natom.and.ipert2>0.and.ipert2<=natom) then
   if (nspinor==1) then
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       if (iatom/=ipert1.or.iatom/=ipert2) cycle ! To move atom "ipert" does not change projectors of other atoms
       cplex_rhoij=pawrhoij(iatom)%cplex_rhoij
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         d1cpj0(1:2,1)=cwaveprj0_pert1(iatm,1)%dcp(1:2,1,jlmn)   ! < p_j^(pert1) | Psi^(0) >
         d2cpj0(1:2,1)=cwaveprj0_pert2(iatm,1)%dcp(1:2,1,jlmn)   ! < p_j^(pert2) | Psi^(0) >
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           d1cpi0(1:2,1)=cwaveprj0_pert1(iatm,1)%dcp(1:2,1,ilmn) ! < p_i^(pert1) | Psi^(0) >
           d2cpi0(1:2,1)=cwaveprj0_pert2(iatm,1)%dcp(1:2,1,ilmn) ! < p_i^(pert2) | Psi^(0) >
           ro11_re=zero
           do iplex=1,cplex
             ro11_re=ro11_re+d1cpi0(iplex,1)*d2cpj0(iplex,1)+d2cpi0(iplex,1)*d1cpj0(iplex,1)
           end do
           pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol)+weight*ro11_re
           if (compute_impart_cplex) then
             klmn_im=klmn_re+1
             ro11_im=        d1cpi0(1,1)*d2cpj0(2,1)-d1cpi0(2,1)*d2cpj0(1,1)
             ro11_im=ro11_im+d2cpi0(1,1)*d1cpj0(2,1)-d2cpi0(2,1)*d1cpj0(1,1)
             pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol)+weight*ro11_im
           end if
         end do
       end do
     end do
   else ! nspinor=2
     message="paw_dfptnl_accrhoij is not implemented for nspinor=2!"
     MSG_BUG(message)
   end if
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine paw_dfptnl_accrhoij
!!***

!----------------------------------------------------------------------

END MODULE m_paw_dfptnl
!!***

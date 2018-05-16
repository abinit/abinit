!{\src2tex{textfont=tt}}
!!****f* ABINIT/paw_dfptnl_energy
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

!! COPYRIGHT
!! Copyright (C) 2018-2018 ABINIT group (LB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
 use m_pawrhoij,   only : pawrhoij_type
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
 integer,intent(in) :: pawprtvol,pawxcdev
 integer,optional,intent(in) :: comm_atom
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(out) :: d3exc(2)
 type(paw_an_type),intent(in) :: paw_an0(my_natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij_1(natom)
 type(pawrhoij_type),intent(in) :: pawrhoij_2(natom)
 type(pawrhoij_type),intent(in) :: pawrhoij_3(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: cplex_1,cplex_2,cplex_3,iatom,iatom_tot,itypat
 integer :: lm_size_all,mesh_size,my_comm_atom,npts,nspden,nzlmopt
 integer :: opt_compch,usecore,usetcore,usexcnhat
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: compch,d3exc1_iat(2)
! character(len=500) :: msg
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
   MSG_BUG("paw_dfptnl_energy is not implemented for pawxcdev /=0")
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!!Various inits
 opt_compch=0;!optvxc=1;optexc=3
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
   cplex_1=pawrhoij_1(iatom)%cplex
   cplex_2=pawrhoij_2(iatom)%cplex
   cplex_3=pawrhoij_3(iatom)%cplex
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
&                 lmselect_1,lmselect_2,lmselect_3,nhat1_1,nhat1_2,nhat1_3,&
&                 paw_an0(iatom)%nk3xc1,mesh_size,nspden,pawang,pawrad(itypat),&
                  rho1_1,rho1_2,rho1_3,0)
   d3exc = d3exc + d3exc1_iat

   call paw_dfptnl_xc(cplex_1,cplex_2,cplex_3,d3exc1_iat,ixc,paw_an0(iatom)%k3xct1,lm_size_all,&
&                 lmselect_1,lmselect_2,lmselect_3,nhat1_1,nhat1_2,nhat1_3,&
&                 paw_an0(iatom)%nk3xc1,mesh_size,nspden,pawang,pawrad(itypat),&
                  trho1_1,trho1_2,trho1_3,usexcnhat)
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

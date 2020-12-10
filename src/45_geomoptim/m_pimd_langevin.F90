!!****m* ABINIT/m_pimd_langevin
!! NAME
!! m_pimd_langevin
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2011-2020 ABINIT group (GG,MT)
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

module m_pimd_langevin

 use defs_basis
 use m_abicore
 use m_errors
 use m_pimd
 use m_random_zbq

 use m_symtk,     only : matr3inv
 use m_geometry,  only : xcart2xred, xred2xcart


 implicit none

 private
!!***

 public :: pimd_langevin_npt
 public :: pimd_langevin_nvt
!!***

contains
!!***

!!****f* ABINIT/pimd_langevin_npt
!! NAME
!! pimd_langevin_npt
!!
!! FUNCTION
!! Predicts new positions in Path Integral Molecular Dynamics using Langevin thermostat in the NPT ensemble.
!! Given the positions at time t and t-dtion, an estimation of the velocities at time t,
!! the forces and an estimation of the stress at time t, and an estimation of the cell at time t,
!! computes in the Path Integral Molecular Dynamics framework the new positions at time t+dtion,
!! computes self-consistently the velocities, the stress and the cell at time t and produces
!! an estimation of the velocities, stress and new cell at time t+dtion
!!
!! INPUTS
!!  etotal(trotter)=electronic total energy for all images
!!  itimimage=number of the current time for image propagation (itimimage+1 is to be predicted here)
!!  natom=dimension of vel_timimage and xred_timimage
!!  pimd_param=datastructure that contains all the parameters necessary to Path-Integral MD
!!  prtvolimg=printing volume
!!  rprimd(3,3)=dimensionless unit cell vectors (common to all images) at time t (present time step)
!!  rprimd_prev(3,3)=dimensionless unit cell vectors (common to all images) at time t-dt (previous time step)
!!  stressin(3,3,trotter)=electronic stress tensor for each image
!!  trotter=Trotter number (total number of images)
!!  volume=voume of unit cell (common to all images)
!!  xred(3,natom,trotter)=reduced coordinates of atoms for all images at time t (present time step)
!!  xred_prev(3,natom,trotter)=reduced coordinates of atoms for all images at time t-dt (previous time step)
!!
!! OUTPUT
!!  rprimd_next(3,3)=dimensionless unit cell vectors (common to all images) at time t+dt (next time step)
!!  xred_next(3,natom,trotter)=reduced coordinates of atoms for all images at time t+dt (next time step)
!!
!! SIDE EFFECTS
!!  forces(3,natom,trotter)=forces over atoms for all images
!!    at input,  electronic forces
!!    at output, electronic forces + quantum spring contribution
!!  vel(3,natom,trotter)=velocies of atoms for all images
!!    at input,  values at time t
!!    at output, values at time t+dt
!!  vel_cell(3,3,trotter)=time derivative of cell parameters
!!    at input,  values at time t
!!    at output, values at time t+dt
!!
!! NOTES
!!  Here follows PIMD in the NPT ensemble within the Langevin barostat algorithm
!!  of Quigley and Probert: J. Chem. Phys. 120, 11432 (2004) [[cite:Quigley2004]]
!!  and Comput. Phys. Comm. 169, 322 (2005) [[cite:Quigley2005]]
!!
!! PARENTS
!!      m_predict_pimd
!!
!! CHILDREN
!!      pimd_apply_constraint,pimd_coord_transform,pimd_energies
!!      pimd_force_transform,pimd_forces,pimd_initvel,pimd_langevin_forces
!!      pimd_langevin_random,pimd_langevin_random_init,pimd_langevin_random_qtb
!!      pimd_mass_spring,pimd_predict_taylor,pimd_predict_verlet,pimd_print
!!      pimd_stresses,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine pimd_langevin_npt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&
&          rprimd,rprimd_next,rprimd_prev,stressin,trotter,vel,vel_cell,&
&          volume,xred,xred_next,xred_prev)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimimage,natom,prtvolimg,trotter
 real(dp),intent(in) :: volume
 type(pimd_type),intent(in) :: pimd_param
!arrays
 real(dp),intent(in) :: etotal(trotter),rprimd(3,3),rprimd_prev(3,3),stressin(3,3,trotter)
 real(dp),intent(in),target :: xred(3,natom,trotter),xred_prev(3,natom,trotter)
 real(dp),intent(out) :: rprimd_next(3,3),xred_next(3,natom,trotter)
 real(dp),intent(inout) :: forces(3,natom,trotter),vel(3,natom,trotter),vel_cell(3,3,trotter)

!Local variables-------------------------------
!Options
!        Option for the Langevin algorithm correction
 integer,parameter :: ilangevin=0
!        The following option forces the total of forces to be zero
!         It prevents the translation of the center of mass
!         If it is zero, no constraint on mass center is applied
 integer,parameter :: zeroforce=1
!        Tolerance for the SC cycle
 real(dp),parameter :: tolerance=tol7

!scalars
 integer :: idum=-5
 integer :: constraint,iatom,ii,iimage,irestart,jj,ndof,prtstress
 real(dp) :: dtion,eharm,eharm2,epot,friction,frictionbar,initemp,kt,rescale_temp,scalebar
 real(dp) :: temperature1,temperature2,temp2_prev,thermtemp,tol,tracepg,wg
!arrays
 real, parameter :: identity(3,3)=reshape((/(one,(zero,ii=1,3),jj=1,2),one/),(/3,3/))
 real(dp) :: aleabar(3,3),constraint_output(2),ddh(3,3),diffstress(3,3)
 real(dp) :: dstrhh(3,3),fg(3,3),invrprimd(3,3)
 real(dp) :: langev_bar(3,3),pg(3,3),pgdh(3,3),stress_pimd(3,3,3),strtarget(6),tmp(3,3)
 real(dp),allocatable :: alea(:,:,:),forces_orig(:,:,:),forces_pimd(:,:,:),forces_pimd_red(:,:)
 real(dp),allocatable :: fsup(:,:),hxredpoint(:,:,:),inertmass(:),langev(:,:)
 real(dp),allocatable ::  mass(:,:),quantummass(:),spring(:,:)
 real(dp),allocatable :: xcart(:,:,:),xcart_next(:,:,:),xcart_prev(:,:,:)
 real(dp),allocatable :: xredpoint(:,:,:)

! *************************************************************************

!############# Initializations ###########################

!Allocation of local arrays
 ABI_ALLOCATE(xcart,(3,natom,trotter))
 ABI_ALLOCATE(xcart_prev,(3,natom,trotter))
 ABI_ALLOCATE(xcart_next,(3,natom,trotter))
 ABI_ALLOCATE(forces_orig,(3,natom,trotter))
 ABI_ALLOCATE(forces_pimd,(3,natom,trotter))
 ABI_ALLOCATE(inertmass,(natom))
 ABI_ALLOCATE(quantummass,(natom))

!Fill in the local variables
 ndof=3*natom*trotter
 rescale_temp=one;if(zeroforce==1)rescale_temp=dble(ndof)/dble(ndof-3)
 quantummass(1:natom)=pimd_param%amu   (pimd_param%typat(1:natom))*amu_emass
 inertmass  (1:natom)=pimd_param%pimass(pimd_param%typat(1:natom))*amu_emass
 initemp=pimd_param%mdtemp(1)/rescale_temp; thermtemp=pimd_param%mdtemp(2)
 dtion=pimd_param%dtion
 kt=thermtemp*kb_HaK
 friction=pimd_param%vis
 wg=pimd_param%bmass
 strtarget(:)=pimd_param%strtarget(:) ! imposed stress tensor
 frictionbar=pimd_param%friction      ! friction coeff of barostat
 scalebar=sqrt(two*frictionbar*wg*kt/dtion)
 forces_orig=forces

!Masses and spring constants
 ABI_ALLOCATE(mass,(natom,1))
 ABI_ALLOCATE(spring,(natom,1))
 call pimd_mass_spring(inertmass,kt,mass,natom,quantummass,spring,0,trotter)

!Initialize random forces
 ABI_ALLOCATE(alea,(3,natom,trotter))
 ABI_ALLOCATE(langev,(natom,trotter))
 langev(:,1)=sqrt(two*friction*inertmass(:)*kt/dtion)
 if(ilangevin==1)then
   langev(:,1)=langev(:,1)*sqrt(one-(friction*dtion/(two*inertmass(:))))
 end if

!Random number generator initialization
 if(itimimage<=1) then
   call pimd_langevin_random_init(pimd_param%irandom,idum)
 end if

!Compute cartesian coordinates
 do iimage=1,trotter
   call xred2xcart(natom,rprimd,xcart     (:,:,iimage),xred(:,:,iimage))
   call xred2xcart(natom,rprimd,xcart_prev(:,:,iimage),xred_prev(:,:,iimage))
 end do

!Determine if it is a restart or not
!If this is a calculation from scratch,generate random distribution of velocities
 irestart=1;if (itimimage==1) irestart=pimd_is_restart(mass,vel,vel_cell)

!Initialize derivatives
 if (mod(irestart,10)==0) then
   call pimd_initvel(idum,mass,natom,initemp,trotter,vel,pimd_param%constraint,pimd_param%wtatcon)
 end if
!vel_cell does not depend on Trotter...
 ddh=vel_cell(:,:,1);if (irestart<10) ddh=zero

 if (itimimage<=1) then

!  ========================= FIRST TIME STEP =======================================

   ABI_ALLOCATE(hxredpoint,(3,natom,trotter))
   ABI_ALLOCATE(xredpoint,(3,natom,trotter))
   ABI_ALLOCATE(forces_pimd_red,(3,natom))
   ABI_ALLOCATE(fsup,(3,natom))

   do iimage=1,trotter
     hxredpoint(:,:,iimage)=vel(:,:,iimage) - matmul(ddh(:,:),xred(:,:,iimage))
   end do
   call matr3inv(rprimd,invrprimd)
   do iimage=1,trotter
     xredpoint(:,:,iimage)=matmul(invrprimd(:,:),hxredpoint(:,:,iimage))
   end do

!  Generate random numbers
   call pimd_langevin_random(alea,pimd_param%irandom,idum,langev,mass,natom,trotter,zeroforce)

!  Compute PIMD and Langevin contributions to forces
   call pimd_forces(forces,natom,spring,0,trotter,xcart)
   call pimd_langevin_forces(alea,forces,forces_pimd,friction,langev,mass,natom,trotter,hxredpoint)
   call pimd_apply_constraint(pimd_param%constraint,constraint_output,forces_pimd,&
&   mass,natom,trotter,pimd_param%wtatcon,xcart)
   tmp=matmul(invrprimd,ddh)
   pg=wg*matmul(ddh,invrprimd)
   tracepg=pg(1,1)+pg(2,2)+pg(3,3)

!  Taylor algorithm
   do iimage=1,trotter
     call xcart2xred(natom,rprimd,forces_pimd(:,:,iimage),forces_pimd_red)
     fsup(:,:)=matmul(tmp,xredpoint(:,:,iimage))
     do iatom=1,natom
       xred_next(:,iatom,iimage)=xred(:,iatom,iimage)+dtion*xredpoint(:,iatom,iimage) + &
&       half*(  &
&       forces_pimd_red(:,iatom)-two*inertmass(iatom)*fsup(:,iatom)- &
&       (tracepg*inertmass(iatom)*xredpoint(:,iatom,iimage)/(wg*dble(ndof))) &
&       )*dtion*dtion/inertmass(iatom)
     end do
   end do

!  predict rprimd at time t+dt from taylor algorithm
   call pimd_langevin_random_bar(aleabar,pimd_param%irandom,idum)
   langev_bar=matmul(aleabar,rprimd)*scalebar
   call pimd_stresses(mass,natom,quantummass,stress_pimd,stressin,thermtemp,thermtemp,trotter,hxredpoint,volume,xcart)
   diffstress=pimd_diff_stress(stress_pimd,strtarget)
   dstrhh=matmul(diffstress,rprimd)
   pgdh=matmul(pg,ddh)
   temperature1=pimd_temperature(mass,hxredpoint)*rescale_temp
   fg(:,:)=volume*dstrhh(:,:)+pgdh(:,:)+temperature1*kb_HaK*rprimd(:,:)

   rprimd_next(:,:)=rprimd(:,:) + dtion*ddh(:,:) + half*(  &
&   fg(:,:)-wg*frictionbar*ddh(:,:)+langev_bar  &
&   )*dtion*dtion/wg

!  Recompute xcart_next
   do iimage=1,trotter
     call xred2xcart(natom,rprimd_next,xcart_next(:,:,iimage),xred_next(:,:,iimage))
   end do

!  Compute stress tensor at t from virial theorem
   call pimd_stresses(mass,natom,quantummass,stress_pimd,stressin,thermtemp,thermtemp,trotter,vel,volume,xcart)

!  Translate from pressure to stress for print
   stress_pimd=-stress_pimd

!  Compute temperature at current step
   temperature1=pimd_temperature(mass,vel)*rescale_temp

!  Estimate the velocities at t+dt/2
   vel=(xcart_next-xcart)/dtion
   ddh=(rprimd_next-rprimd)/dtion

!  Compute new temperature
   temperature2=pimd_temperature(mass,vel)*rescale_temp

   ABI_DEALLOCATE(hxredpoint)
   ABI_DEALLOCATE(xredpoint)
   ABI_DEALLOCATE(forces_pimd_red)
   ABI_DEALLOCATE(fsup)

 else

!  ========================= OTHER TIME STEPS ======================================

!  Additional allocations
   ABI_ALLOCATE(hxredpoint,(3,natom,trotter))
   ABI_ALLOCATE(xredpoint,(3,natom,trotter))
   ABI_ALLOCATE(forces_pimd_red,(3,natom))
   ABI_ALLOCATE(fsup,(3,natom))

!  first estimation of ddh, pg and its trace:
   call matr3inv(rprimd,invrprimd)
   pg=wg*matmul(ddh,invrprimd)
   tracepg=pg(1,1)+pg(2,2)+pg(3,3)

!  Momenta hxredpoint = H ds/dt: estimation
   if (itimimage==2) then
     hxredpoint=vel
   else
     do iimage=1,trotter
       hxredpoint(:,:,iimage)=matmul(rprimd,vel(:,:,iimage))
!      because vel in entrance is a scaled velocity (ds/dt)
     end do
   end if

!  Compute temperature at t
   temperature1=pimd_temperature(mass,hxredpoint)*rescale_temp

!  Generate random numbers
   call pimd_langevin_random(alea,pimd_param%irandom,idum,langev,mass,natom,trotter,zeroforce)

!  Generate random numbers for the barostat
   call pimd_langevin_random_bar(aleabar,pimd_param%irandom,idum)
   langev_bar=matmul(aleabar,rprimd)*scalebar

!  Compute PIMD and Langevin contributions to forces
   call pimd_forces(forces,natom,spring,0,trotter,xcart)
   call pimd_langevin_forces(alea,forces,forces_pimd,friction,langev,mass,natom,trotter,hxredpoint)
   call pimd_apply_constraint(pimd_param%constraint,constraint_output,forces_pimd,&
&   mass,natom,trotter,pimd_param%wtatcon,xcart)

!  Compute difference between instantaneous stress and imposed stress (barostat)
   call pimd_stresses(mass,natom,quantummass,stress_pimd,stressin,thermtemp,thermtemp,trotter,hxredpoint,volume,xcart)

   diffstress=pimd_diff_stress(stress_pimd,strtarget)

!  Compute "force" on supercell vectors
   dstrhh=matmul(diffstress,rprimd)
   pgdh=matmul(pg,ddh)
   fg(:,:)=volume*dstrhh(:,:)+pgdh(:,:)+temperature1*kb_HaK*rprimd(:,:)

!  Evolve the supercell (fist estimation)
   rprimd_next=two*rprimd-rprimd_prev+(fg-wg*frictionbar*ddh+langev_bar)*dtion*dtion/wg

!  Evolve atomic positions (first estimation)
   tmp=matmul(invrprimd,ddh)
   do iimage=1,trotter
     call xcart2xred(natom,rprimd,hxredpoint(:,:,iimage),xredpoint(:,:,iimage))
     fsup(:,:)=matmul(tmp,xredpoint(:,:,iimage))
     call xcart2xred(natom,rprimd,forces_pimd(:,:,iimage),forces_pimd_red)

     do iatom=1,natom
       xred_next(:,iatom,iimage)= &
&       two*xred(:,iatom,iimage) - xred_prev(:,iatom,iimage) &
&       +(forces_pimd_red(:,iatom)-two*inertmass(iatom)*fsup(:,iatom) &
&       -tracepg*inertmass(iatom)*xredpoint(:,iatom,iimage)/(wg*dble(ndof))) &
&       *dtion*dtion/inertmass(iatom)
     end do
   end do

!  Self-consistent loop
   temperature2=pimd_temperature(mass,xredpoint)*rescale_temp
   temp2_prev=temperature2; tol=one

   do while (tol>tolerance)
!    Reestimate dH/dt at t
     ddh(:,:)=(rprimd_next(:,:)-rprimd_prev(:,:))/(two*dtion)

!    Reestimate the scaled velocities at t
     do iimage=1,trotter
       xredpoint(:,:,iimage)=(xred_next(:,:,iimage)-xred_prev(:,:,iimage))/(two*dtion)
       call xred2xcart(natom,rprimd,hxredpoint(:,:,iimage),xredpoint(:,:,iimage))
     end do
!    Reestimate the forces
     call pimd_langevin_forces(alea,forces,forces_pimd,friction,langev,mass,natom,trotter,hxredpoint)
     call pimd_apply_constraint(pimd_param%constraint,constraint_output,forces_pimd,&
&     mass,natom,trotter,pimd_param%wtatcon,xcart)
!    Compute variation of temperature (to check convergence of SC loop)
     temperature2=pimd_temperature(mass,xredpoint)*rescale_temp
     tol=dabs(temperature2-temp2_prev)/dabs(temp2_prev)
     temp2_prev=temperature2
!    Recompute the temperature
     temperature2=pimd_temperature(mass,hxredpoint)*rescale_temp
!    Recompute pg
     pg=wg*matmul(ddh,invrprimd)
     tracepg=pg(1,1)+pg(2,2)+pg(3,3)
!    Recompute difference between instantaneous stress and imposed stress (barostat)
     call pimd_stresses(mass,natom,quantummass,stress_pimd,stressin,thermtemp,thermtemp,trotter,hxredpoint,volume,xcart)

     diffstress=pimd_diff_stress(stress_pimd,strtarget)

!    Recompute "force" on supercell vectors
     dstrhh=matmul(diffstress,rprimd)
     pgdh=matmul(pg,ddh)
     fg(:,:)=volume*diffstress(:,:)+pgdh(:,:)+temperature2*kb_HaK*rprimd(:,:)
!    Evolve the supercell (better estimation)
     rprimd_next=two*rprimd-rprimd_prev+(fg-wg*frictionbar*ddh+langev_bar)*dtion*dtion/wg

!    Evolve atomic positions (better estimation):
     tmp=matmul(invrprimd,ddh)
     do iimage=1,trotter
       call xcart2xred(natom,rprimd,hxredpoint(:,:,iimage),xredpoint(:,:,iimage))
       fsup(:,:)=matmul(tmp,xredpoint(:,:,iimage))
       call xcart2xred(natom,rprimd,forces_pimd(:,:,iimage),forces_pimd_red)
       do iatom=1,natom
         xred_next(:,iatom,iimage)= &
&         two*xred(:,iatom,iimage) - xred_prev(:,iatom,iimage) &
&         +(forces_pimd_red(:,iatom)-two*inertmass(iatom)*fsup(:,iatom) &
&         -tracepg*inertmass(iatom)*xredpoint(:,iatom,iimage)/(wg*dble(ndof))) &
&         *dtion*dtion/inertmass(iatom)
       end do
     end do
   end do ! End self-consistent loop

!  Computation of true temperature from true velocities at t
   do iimage=1,trotter
     vel(:,:,iimage)=hxredpoint(:,:,iimage)+matmul(ddh,xred(:,:,iimage))
   end do
   temperature2=pimd_temperature(mass,vel)*rescale_temp

!  Computation of the real stress tensor at t
   call pimd_stresses(mass,natom,quantummass,stress_pimd,stressin,thermtemp,thermtemp,trotter,vel,volume,xcart)

!  translate from pressure to stress
   stress_pimd=-stress_pimd

!  Deallocations (Verlet algo)
   ABI_DEALLOCATE(xredpoint)
   ABI_DEALLOCATE(hxredpoint)
   ABI_DEALLOCATE(forces_pimd_red)
   ABI_DEALLOCATE(fsup)

 end if ! itimimage==1

!############# Final operations ############################

!Compute contributions to energy
 call pimd_energies(eharm,eharm2,epot,etotal,forces_orig,natom,spring,trotter,xcart)

!Print messages
 prtstress=1
 call pimd_print(constraint,constraint_output,&
& eharm,eharm2,epot,forces_pimd,inertmass,irestart,&
& itimimage,kt,natom,pimd_param%optcell,prtstress,prtvolimg,rprimd,&
& stress_pimd,temperature1,temperature2,&
& pimd_param%traj_unit,trotter,vel,ddh,xcart,xred)

 if (itimimage>1) then
!  Estimation of ds/dt at t+dt
   vel = (three*xred_next - four*xred + xred_prev)/(two * dtion)
   ddh = (three*rprimd_next - four*rprimd + rprimd_prev)/(two * dtion)
 end if

!Come back to reduced coordinates
 if(itimimage<=1) then
   do iimage=1,trotter
     call xcart2xred(natom,rprimd,xcart_next(:,:,iimage),xred_next(:,:,iimage))
   end do
 end if

!Return cell velocities (does not depend on Trotter)
 do iimage=1,trotter
   vel_cell(:,:,iimage)=ddh(:,:)
 end do

!Free memory
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xcart_prev)
 ABI_DEALLOCATE(xcart_next)
 ABI_DEALLOCATE(forces_orig)
 ABI_DEALLOCATE(forces_pimd)
 ABI_DEALLOCATE(inertmass)
 ABI_DEALLOCATE(quantummass)
 ABI_DEALLOCATE(mass)
 ABI_DEALLOCATE(spring)
 ABI_DEALLOCATE(alea)
 ABI_DEALLOCATE(langev)

end subroutine pimd_langevin_npt
!!***

!!****f* ABINIT/pimd_langevin_nvt
!! NAME
!! pimd_langevin_nvt
!!
!! FUNCTION
!! Predicts new positions in Path Integral Molecular Dynamics using Langevin thermostat in the NVT ensemble.
!! Given the positions at time t and t-dtion, an estimation of the velocities at time t,
!! the forces and an estimation of the stress at time t, and an estimation of the cell at time t,
!! computes in the Path Integral Molecular Dynamics framework the new positions at time t+dtion,
!! computes self-consistently the velocities, the stress and the cell at time t and produces
!! an estimation of the velocities, stress and new cell at time t+dtion
!!
!! INPUTS
!!  etotal(trotter)=electronic total energy for all images
!!  itimimage=number of the current time for image propagation (itimimage+1 is to be predicted here)
!!  natom=dimension of vel_timimage and xred_timimage
!!  pimd_param=datastructure that contains all the parameters necessary to Path-Integral MD
!!  prtvolimg=printing volume
!!  rprimd(3,3)=dimensionless unit cell vectors (common to all images)
!!  stressin(3,3,trotter)=electronic stress tensor for each image
!!  trotter=Trotter number (total number of images)
!!  volume=voume of unit cell (common to all images)
!!  xred(3,natom,trotter)=reduced coordinates of atoms for all images at time t (present time step)
!!  xred_prev(3,natom,trotter)=reduced coordinates of atoms for all images at time t-dt (previous time step)
!!
!! OUTPUT
!!  xred_next(3,natom,trotter)=reduced coordinates of atoms for all images at time t+dt (next time step)
!!
!! SIDE EFFECTS
!!  forces(3,natom,trotter)=forces over atoms for all images
!!    at input,  electronic forces
!!    at output, electronic forces + quantum spring contribution
!!  vel(3,natom,trotter)=velocies of atoms for all images
!!    at input,  values at time t
!!    at output, values at time t+dt
!!
!! NOTES
!!   See Quigley,Probert, JCP 120, 11432 (2004) [[cite:Quigley2004]], part III
!!
!! PARENTS
!!      m_predict_pimd
!!
!! CHILDREN
!!      pimd_apply_constraint,pimd_coord_transform,pimd_energies
!!      pimd_force_transform,pimd_forces,pimd_initvel,pimd_langevin_forces
!!      pimd_langevin_random,pimd_langevin_random_init,pimd_langevin_random_qtb
!!      pimd_mass_spring,pimd_predict_taylor,pimd_predict_verlet,pimd_print
!!      pimd_stresses,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine pimd_langevin_nvt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&
&                            rprimd,stressin,trotter,vel,volume,xred,xred_next,xred_prev)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimimage,natom,prtvolimg,trotter
 real(dp),intent(in) :: volume
 type(pimd_type),intent(in) :: pimd_param
!arrays
 real(dp),intent(in) :: etotal(trotter),rprimd(3,3),stressin(3,3,trotter)
 real(dp),intent(in),target :: xred(3,natom,trotter),xred_prev(3,natom,trotter)
 real(dp),intent(out) :: xred_next(3,natom,trotter)
 real(dp),intent(inout) :: forces(3,natom,trotter),vel(3,natom,trotter)

!Local variables-------------------------------
!Options
!        Option for the Langevin algorithm correction
 integer,parameter :: ilangevin=0
!        Tolerance for the SC cycle
 real(dp),parameter :: tolerance=tol9

!scalars
 integer :: idum=-5
 integer :: iimage,irestart,ndof,pitransform,prtstress,use_qtb,zeroforce
 real(dp) :: dtion,eharm,eharm2,epot,friction,initemp,kt,kt_,rescale_temp
 real(dp) :: temperature1,temperature2,temp2_prev,thermtemp,tol
!arrays
 real(dp) :: constraint_output(2),spring_prim(natom),stress_pimd(3,3,3),vel_cell(3,3)
 real(dp),allocatable :: alea(:,:,:),forces_orig(:,:,:),forces_pimd(:,:,:)
 real(dp),allocatable :: inertmass(:),langev(:,:),mass(:,:),quantummass(:),spring(:,:)
 real(dp),allocatable :: xcart(:,:,:),xcart_next(:,:,:),xcart_prev(:,:,:)

! *************************************************************************

 if (pimd_param%use_qtb==1.and.pimd_param%qtb_file_unit<=0) then
   ABI_BUG('piqtb_force not open!')
 end if

!############# Initializations ###########################

 pitransform=pimd_param%pitransform

!The following option forces the total of forces to be zero
!It prevents the translation of the center of mass
!If it is zero, no constraint on mass center is applied
 zeroforce=1
 if(pitransform==1) zeroforce=0
 if(pitransform==2) zeroforce=0
 if(pimd_param%constraint==1) zeroforce=0

!Allocation of local arrays
 ABI_ALLOCATE(xcart,(3,natom,trotter))
 ABI_ALLOCATE(xcart_prev,(3,natom,trotter))
 ABI_ALLOCATE(xcart_next,(3,natom,trotter))
 ABI_ALLOCATE(forces_orig,(3,natom,trotter))
 ABI_ALLOCATE(forces_pimd,(3,natom,trotter))
 ABI_ALLOCATE(inertmass,(natom))
 ABI_ALLOCATE(quantummass,(natom))

!Fill in the local variables
 use_qtb=pimd_param%use_qtb
 ndof=3*natom*trotter
 rescale_temp=one;if(zeroforce==1)rescale_temp=dble(ndof)/dble(ndof-3)
 quantummass(1:natom)=pimd_param%amu   (pimd_param%typat(1:natom))*amu_emass
 inertmass  (1:natom)=pimd_param%pimass(pimd_param%typat(1:natom))*amu_emass
 if(pitransform==1) inertmass=quantummass !compulsory for good definition of normal mode masses
 if(pitransform==2) inertmass=quantummass !compulsory for good definition of staging masses
 initemp=pimd_param%mdtemp(1)/rescale_temp;thermtemp=pimd_param%mdtemp(2)
 friction=pimd_param%vis;dtion=pimd_param%dtion
 kt=thermtemp*kb_HaK
 forces_orig=forces

!Masses and spring constants
 select case(pitransform)
 case(0)
   ABI_ALLOCATE(mass,(natom,1))   ! This second dimension is needed
   ABI_ALLOCATE(spring,(natom,1))
   ABI_ALLOCATE(langev,(natom,1))
 case(1,2)
   ABI_ALLOCATE(mass,(natom,trotter))
   ABI_ALLOCATE(spring,(natom,trotter))
   ABI_ALLOCATE(langev,(natom,trotter))
 end select
 spring_prim(:)=quantummass(:)*dble(trotter)*kt*kt
 call pimd_mass_spring(inertmass,kt,mass,natom,quantummass,spring,pitransform,trotter)

!Initialize random forces
 ABI_ALLOCATE(alea,(3,natom,trotter))
 if (use_qtb==0) then
   langev(:,:)=sqrt(two*friction*mass(:,:)*kt/dtion)
 else
   langev(:,:)=sqrt(two*friction*mass(:,:))
 end if

!Random number generator initialization
 if(itimimage<=1) then
   call pimd_langevin_random_init(pimd_param%irandom,idum)
 end if

!Compute cartesian coordinates
 do iimage=1,trotter
   call xred2xcart(natom,rprimd,xcart     (:,:,iimage),xred(:,:,iimage))
   call xred2xcart(natom,rprimd,xcart_prev(:,:,iimage),xred_prev(:,:,iimage))
 end do

!Determine if it is a restart or not
!If this is a calculation from scratch,generate random distribution of velocities
 irestart=1;if (itimimage==1) irestart=pimd_is_restart(mass,vel)
 if (irestart==0) then
   call pimd_initvel(idum,mass,natom,initemp,trotter,vel,pimd_param%constraint,pimd_param%wtatcon)
 end if

!Compute temperature at t
 temperature1=pimd_temperature(mass,vel)*rescale_temp

!################## Images evolution #####################

!Generate random numbers
 if (use_qtb==0) then
   call pimd_langevin_random(alea,pimd_param%irandom,idum,langev,mass,natom,trotter,zeroforce)
 else
   call pimd_langevin_random_qtb(alea,langev,mass,natom,pimd_param%qtb_file_unit,trotter,zeroforce)
 end if

!Compute PIMD and Langevin contributions to forces
 call pimd_coord_transform(xcart,1,natom,pitransform,trotter)
 call pimd_force_transform(forces,1,natom,pitransform,trotter) !compute staging forces
 call pimd_forces(forces,natom,spring,pitransform,trotter,xcart)
 call pimd_langevin_forces(alea,forces,forces_pimd,friction,langev,mass,natom,trotter,vel)
 call pimd_apply_constraint(pimd_param%constraint,constraint_output,forces_pimd,&
& mass,natom,trotter,pimd_param%wtatcon,xcart)

!Compute atomic positions at t+dt
 if (itimimage<=1) then

!  === 1st time step: single Taylor algorithm
!  Predict positions
   call pimd_predict_taylor(dtion,forces_pimd,mass,natom,trotter,&
&   vel,xcart,xcart_next)

!  Estimate the velocities at t+dt/2
   vel=(xcart_next-xcart)/dtion

!  Compute new temperature
   temperature2=pimd_temperature(mass,vel)*rescale_temp

 else

!  === Other time steps: Verlet algorithm + SC cycle
!  Predict positions
   call pimd_coord_transform(xcart_prev,1,natom,pitransform,trotter)
   call pimd_predict_verlet(dtion,forces_pimd,mass,natom,trotter,&
&   xcart,xcart_next,xcart_prev)
!  Self-consistent loop
   temperature2=pimd_temperature(mass,vel)*rescale_temp
   temp2_prev=temperature2; tol=one
   do while (tol>tolerance)
!    Recompute a (better) estimation of the velocity at time step t
     vel = (xcart_next - xcart_prev) / (two*dtion)
     temperature2=pimd_temperature(mass,vel)*rescale_temp
!    Reestimate the force
     call pimd_langevin_forces(alea,forces,forces_pimd,friction,&
&     langev,mass,natom,trotter,vel)
     call pimd_apply_constraint(pimd_param%constraint,constraint_output,forces_pimd,&
&     mass,natom,trotter,pimd_param%wtatcon,xcart)
!    Compute new positions
     call pimd_predict_verlet(dtion,forces_pimd,mass,natom,trotter,&
&     xcart,xcart_next,xcart_prev)

!    Compute variation of temperature (to check convergence of SC loop)
     tol=dabs(temperature2-temp2_prev)/dabs(temp2_prev)
     temp2_prev=temperature2

   end do ! End self-consistent loop

 end if ! itimimage==1

 call pimd_coord_transform(xcart_next,-1,natom,pitransform,trotter)
 call pimd_coord_transform(xcart,-1,natom,pitransform,trotter)
 call pimd_coord_transform(xcart_prev,-1,natom,pitransform,trotter)

!Compute contributions to energy
 call pimd_energies(eharm,eharm2,epot,etotal,forces_orig,natom,spring_prim,trotter,xcart)

!Compute stress tensor at t
 if (use_qtb==0) then
   call pimd_stresses(mass,natom,quantummass,stress_pimd,stressin,thermtemp,thermtemp,trotter,vel,volume,xcart)
 else
   call pimd_stresses(mass,natom,quantummass,stress_pimd,stressin,temperature2,thermtemp,trotter,vel,volume,xcart)
 end if
 stress_pimd=-stress_pimd ! Translate pressure to stress

!############# Final operations ############################

!Print messages
 vel_cell=zero;prtstress=1;if (prtvolimg>=2) prtstress=0
 kt_=kt;if (use_qtb==1) kt_=temperature2*kb_HaK
 call pimd_print(pimd_param%constraint,constraint_output,&
& eharm,eharm2,epot,forces_pimd,inertmass,irestart,&
& itimimage,kt_,natom,pimd_param%optcell,prtstress,prtvolimg,rprimd,&
& stress_pimd,temperature1,temperature2,&
& pimd_param%traj_unit,trotter,vel,vel_cell,xcart,xred)

!If possible, estimate the (transformed) velocities at t+dt
 if (itimimage>1) then
   call pimd_coord_transform(xcart_next,1,natom,pitransform,trotter)
   call pimd_coord_transform(xcart,1,natom,pitransform,trotter)
   call pimd_coord_transform(xcart_prev,1,natom,pitransform,trotter)
   vel = (three*xcart_next - four*xcart + xcart_prev)/(two * dtion)
   call pimd_coord_transform(xcart_next,-1,natom,pitransform,trotter)
   call pimd_coord_transform(xcart,-1,natom,pitransform,trotter)
   call pimd_coord_transform(xcart_prev,-1,natom,pitransform,trotter)
 end if

!Come back to reduced coordinates
 do iimage=1,trotter
   call xcart2xred(natom,rprimd,xcart_next(:,:,iimage),xred_next(:,:,iimage))
 end do

!Free memory
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xcart_prev)
 ABI_DEALLOCATE(xcart_next)
 ABI_DEALLOCATE(forces_orig)
 ABI_DEALLOCATE(forces_pimd)
 ABI_DEALLOCATE(inertmass)
 ABI_DEALLOCATE(quantummass)
 ABI_DEALLOCATE(mass)
 ABI_DEALLOCATE(spring)
 ABI_DEALLOCATE(alea)
 ABI_DEALLOCATE(langev)

end subroutine pimd_langevin_nvt
!!***

end module m_pimd_langevin
!!***

!!****m* ABINIT/m_pimd_nosehoover
!! NAME
!!   m_pimd_nosehoover
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

module m_pimd_nosehoover

 use defs_basis
 use m_pimd
 use m_abicore

 use m_geometry,  only : xcart2xred, xred2xcart

 implicit none

 private
!!***

 public :: pimd_nosehoover_npt
 public :: pimd_nosehoover_nvt
!!***

contains
!!***

!!****f* ABINIT/pimd_nosehoover_npt
!! NAME
!! pimd_nosehoover_npt
!!
!! FUNCTION
!! Predicts new positions in Path Integral Molecular Dynamics using Nose-Hoover in the NPT ensemble.
!! Given the positions at time t and t-dtion, an estimation of the velocities at time t,
!! the forces and an estimation of the stress at time t, and an estimation of the cell at time t,
!! computes in the Path Integral Molecular Dynamics framework the new positions at time t+dtion,
!! computes self-consistently the velocities, the stress and the cell at time t and produces
!! an estimation of the velocities, stress and new cell at time t+dtion
!! No change of acell and rprim at present.
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
!!    at output, electronic forces + Langevin contribution
!!  vel(3,natom,trotter)=velocies of atoms for all images
!!    at input,  values at time t
!!    at output, values at time t+dt
!!  vel_cell(3,3,trotter)=time derivative of cell parameters
!!    at input,  values at time t
!!    at output, values at time t+dt
!!
!! PARENTS
!!      m_predict_pimd
!!
!! CHILDREN
!!      pimd_apply_constraint,pimd_coord_transform,pimd_energies
!!      pimd_force_transform,pimd_forces,pimd_initvel,pimd_mass_spring
!!      pimd_nosehoover_forces,pimd_nosehoover_propagate,pimd_predict_taylor
!!      pimd_predict_verlet,pimd_print,pimd_stresses,wrtout,xcart2xred
!!      xred2xcart
!!
!! SOURCE

subroutine pimd_nosehoover_npt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&
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
 real(dp),parameter :: tolerance=tol7 ! SCF tolerance
!scalars
 integer :: idum=-5
 integer :: constraint,iimage,irestart,ndof,nnos,pitransform,prtstress
 real(dp) :: dtion,eharm,eharm2,epot,initemp,kt,temperature1,temperature2,thermtemp
!arrays
 real(dp) :: constraint_output(2),ddh(3,3),stress_pimd(3,3,3)
 real(dp),allocatable :: dzeta(:,:,:,:),forces_orig(:,:,:),forces_pimd(:,:,:)
 real(dp),allocatable :: inertmass(:),masseff(:,:),qmass(:),quantummass(:),springeff(:,:)
 real(dp),allocatable :: xcart(:,:,:),xcart_next(:,:,:),xcart_prev(:,:,:),zeta(:)

! *************************************************************************

!############# Initializations ###########################

!Allocation of local arrays
 ABI_MALLOC(xcart,(3,natom,trotter))
 ABI_MALLOC(xcart_prev,(3,natom,trotter))
 ABI_MALLOC(xcart_next,(3,natom,trotter))
 ABI_MALLOC(forces_orig,(3,natom,trotter))
 ABI_MALLOC(forces_pimd,(3,natom,trotter))
 ABI_MALLOC(inertmass,(natom))
 ABI_MALLOC(quantummass,(natom))

!Fill in the local variables
 ndof=3*natom*trotter
 quantummass(1:natom)=pimd_param%amu   (pimd_param%typat(1:natom))*amu_emass
 inertmass  (1:natom)=pimd_param%pimass(pimd_param%typat(1:natom))*amu_emass
 initemp=pimd_param%mdtemp(1);thermtemp=pimd_param%mdtemp(2)
 dtion=pimd_param%dtion;pitransform=pimd_param%pitransform
 kt=thermtemp*kb_HaK
 forces_orig=forces

!Allocation/initialization of local variables used for Nose-Hoover chains
!Associated variables:
!nnos = number of thermostats
!dzeta(3,natom,trotter,nnos) = variables of thermostats, in (atomic time unit)^(-1)
!qmass(nnos) = masses of thermostats
!specific to PIMD: pitransform = coordinate transformation (0:no; 1:normal mode; 2:staging)
 nnos=pimd_param%nnos
 ABI_MALLOC(qmass,(nnos))
 ABI_MALLOC(zeta,(nnos))
 ABI_MALLOC(dzeta,(3,natom,trotter,nnos))
 qmass(1:nnos)=pimd_param%qmass(1:nnos)
 zeta=zero;dzeta=zero

!Compute cartesian coordinates
 do iimage=1,trotter
   call xred2xcart(natom,rprimd,xcart     (:,:,iimage),xred(:,:,iimage))
   call xred2xcart(natom,rprimd,xcart_prev(:,:,iimage),xred_prev(:,:,iimage))
 end do

!Determine if it is a restart or not
!If this is a calculation from scratch,generate random distribution of velocities
 irestart=1;if (itimimage==1) irestart=pimd_is_restart(masseff,vel,vel_cell)

!Initialize derivatives
 if (mod(irestart,10)==0) then
   call pimd_initvel(idum,masseff,natom,initemp,trotter,vel,pimd_param%constraint,pimd_param%wtatcon)
 end if
!vel_cell does not depend on Trotter...
 ddh=vel_cell(:,:,1);if (irestart<10) ddh=zero

!Compute temperature at t
 temperature1=pimd_temperature(masseff,vel)

!################## Images evolution #####################

!This is temporary
 xcart_next=zero
 rprimd_next=rprimd_prev
 temperature2=pimd_temperature(masseff,vel)

!Compute contributions to energy
 call pimd_energies(eharm,eharm2,epot,etotal,forces_orig,natom,springeff,trotter,xcart)

!Compute stress tensor at t from virial theorem
 call pimd_stresses(masseff,natom,quantummass,stress_pimd,stressin,thermtemp,thermtemp,trotter,vel,volume,xcart)


!############# Final operations ############################

!Print messages
 prtstress=1
 call pimd_print(constraint,constraint_output,&
& eharm,eharm2,epot,forces_pimd,inertmass,irestart,&
& itimimage,kt,natom,pimd_param%optcell,prtstress,prtvolimg,rprimd,&
& stress_pimd,temperature1,temperature2,&
& pimd_param%traj_unit,trotter,vel,ddh,xcart,xred)

!Come back to reduced coordinates
 do iimage=1,trotter
   call xcart2xred(natom,rprimd,xcart_next(:,:,iimage),xred_next(:,:,iimage))
 end do

!Return cell velocities (does not depend on Trotter)
 do iimage=1,trotter
   vel_cell(:,:,iimage)=ddh(:,:)
 end do

!Free memory
 ABI_FREE(xcart)
 ABI_FREE(xcart_prev)
 ABI_FREE(xcart_next)
 ABI_FREE(forces_orig)
 ABI_FREE(forces_pimd)
 ABI_FREE(inertmass)
 ABI_FREE(quantummass)
 ABI_FREE(masseff)
 ABI_FREE(springeff)
 ABI_FREE(qmass)
 ABI_FREE(dzeta)
 ABI_FREE(zeta)

end subroutine pimd_nosehoover_npt
!!***

!!****f* ABINIT/pimd_nosehoover_nvt
!! NAME
!! pimd_nosehoover_nvt
!!
!! FUNCTION
!! Predicts new positions in Path Integral Molecular Dynamics using Nose-Hoover in the NVT ensemble.
!! Given the positions at time t and t-dtion, an estimation of the velocities at time t,
!! the forces and an estimation of the stress at time t, and an estimation of the cell at time t,
!! computes in the Path Integral Molecular Dynamics framework the new positions at time t+dtion,
!! computes self-consistently the velocities, the stress and the cell at time t and produces
!! an estimation of the velocities, stress and new cell at time t+dtion
!! No change of acell and rprim at present.
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
!!  Thermization by Nose-Hoover chains according to
!!  Martyna, Klein, Tuckerman, J. Chem. Phys. 97, 2635 (1992) [[cite:Martyna1992]]
!!  Tuckerman, Marx, Klein, Parrinello, J. Chem. Phys. 104, 5579 (1996) [[cite:Tuckerman1996]]
!!
!! PARENTS
!!      m_predict_pimd
!!
!! CHILDREN
!!      pimd_apply_constraint,pimd_coord_transform,pimd_energies
!!      pimd_force_transform,pimd_forces,pimd_initvel,pimd_mass_spring
!!      pimd_nosehoover_forces,pimd_nosehoover_propagate,pimd_predict_taylor
!!      pimd_predict_verlet,pimd_print,pimd_stresses,wrtout,xcart2xred
!!      xred2xcart
!!
!! SOURCE

subroutine pimd_nosehoover_nvt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&
&                              rprimd,stressin,trotter,vel,volume,xred,xred_next,xred_prev)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimimage,natom,prtvolimg,trotter
 real(dp),intent(in) :: volume
 type(pimd_type),intent(inout) :: pimd_param
!arrays
 real(dp),intent(in) :: etotal(trotter),rprimd(3,3),stressin(3,3,trotter)
 real(dp),intent(in),target :: xred(3,natom,trotter),xred_prev(3,natom,trotter)
 real(dp),intent(out) :: xred_next(3,natom,trotter)
 real(dp),intent(inout) :: forces(3,natom,trotter),vel(3,natom,trotter)

!Local variables-------------------------------
!Options
 real(dp),parameter :: tolerance=tol9 ! SCF tolerance
!scalars
 integer :: idum=-5
 integer :: iimage,irestart,ndof,nnos,pitransform,prtstress
 real(dp) :: dtion,eharm,eharm2,epot,initemp,kt
 real(dp) :: temperature1,temperature2,temp2_prev,thermtemp,tol
 character(len=500) :: msg
!arrays
 real(dp) :: constraint_output(2),spring_prim(natom),stress_pimd(3,3,3),vel_cell(3,3)
 real(dp),allocatable :: forces_orig(:,:,:),forces_pimd(:,:,:)
 real(dp),allocatable :: inertmass(:),mass(:,:),qmass(:),quantummass(:),spring(:,:)
 real(dp),allocatable :: xcart(:,:,:),xcart_next(:,:,:),xcart_prev(:,:,:)
 real(dp),allocatable :: dzeta(:,:,:,:),zeta_prev(:,:,:,:),zeta(:,:,:,:)
 real(dp),allocatable :: zeta_next(:,:,:,:)

! *************************************************************************

!############# Initializations ###########################

!Allocation of local arrays
 ABI_MALLOC(xcart,(3,natom,trotter))
 ABI_MALLOC(xcart_prev,(3,natom,trotter))
 ABI_MALLOC(xcart_next,(3,natom,trotter))
 ABI_MALLOC(forces_orig,(3,natom,trotter))
 ABI_MALLOC(forces_pimd,(3,natom,trotter))
 ABI_MALLOC(inertmass,(natom))
 ABI_MALLOC(quantummass,(natom))

!Fill in the local variables
 ndof=3*natom*trotter
 quantummass(1:natom)=pimd_param%amu   (pimd_param%typat(1:natom))*amu_emass
 inertmass  (1:natom)=pimd_param%pimass(pimd_param%typat(1:natom))*amu_emass
 if(pitransform==1) inertmass=quantummass !compulsory for good definition of normal mode masses
 if(pitransform==2) inertmass=quantummass !compulsory for good definition of staging masses
 initemp=pimd_param%mdtemp(1);thermtemp=pimd_param%mdtemp(2)
 dtion=pimd_param%dtion;pitransform=pimd_param%pitransform
 kt=thermtemp*kb_HaK
 forces_orig=forces

!Allocation/initialization of local variables used for Nose-Hoover chains
!Associated variables:
!nnos = number of thermostats
!zeta,zeta_next,zeta_prev(3,natom,trotter,nnos) = variables of thermostats, dzeta, its time derivative
!qmass(nnos) = masses of thermostats
!specific to PIMD: pitransform = coordinate transformation (0:no; 1:normal mode; 2:staging)
 nnos=pimd_param%nnos
 ABI_MALLOC(qmass,(nnos))
 ABI_MALLOC(zeta_prev,(3,natom,trotter,nnos))
 ABI_MALLOC(zeta,(3,natom,trotter,nnos))
 ABI_MALLOC(zeta_next,(3,natom,trotter,nnos))
 ABI_MALLOC(dzeta,(3,natom,trotter,nnos))
!initialization
 qmass(1:nnos)=pimd_param%qmass(1:nnos)
 zeta_prev(:,:,:,:)=pimd_param%zeta_prev(:,:,:,:)
 zeta(:,:,:,:)     =pimd_param%zeta(:,:,:,:)
 dzeta(:,:,:,:)    =pimd_param%dzeta(:,:,:,:) !unuseful to initialize zeta_next

!Masses and spring constants (according to pitransform)
 select case(pitransform)
 case(0)
   ABI_MALLOC(mass,(natom,1))
   ABI_MALLOC(spring,(natom,1))
 case(1,2)
   ABI_MALLOC(mass,(natom,trotter))
   ABI_MALLOC(spring,(natom,trotter))
 end select
 spring_prim(:)=quantummass(:)*dble(trotter)*kt*kt

 call pimd_mass_spring(inertmass,kt,mass,natom,quantummass,spring,pitransform,trotter)

!Recommended value of Nose mass
 write(msg,'(2a,f9.2,3a)') ch10,&
& ' Recommended value of Nose mass is',one/(dble(trotter)*kt),' (atomic units)',ch10,&
& '(see Tuckerman et al, J. Chem. Phys. 104, 5579 (1996))' ! [[cite:Tuckerman1996]]
 call wrtout(std_out,msg,'COLL')

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
 temperature1=pimd_temperature(mass,vel)

!################## Images evolution #####################

!Transform the coordinates and forces (according to pitransform)
 call pimd_coord_transform(xcart,1,natom,pitransform,trotter)
 call pimd_force_transform(forces,1,natom,pitransform,trotter) !compute staging forces
 call pimd_forces(forces,natom,spring,pitransform,trotter,xcart)
 call pimd_nosehoover_forces(dzeta,forces,forces_pimd,mass,natom,nnos,trotter,vel)
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
   temperature2=pimd_temperature(mass,vel)

!  Propagate the thermostat variables
   call pimd_nosehoover_propagate(dtion,dzeta,mass,natom,nnos,qmass,&
&   thermtemp,trotter,vel,zeta,zeta_next,zeta_prev,itimimage,pitransform)

   dzeta=(zeta_next-zeta)/dtion

 else

!  === Other time steps: Verlet algorithm + SC cycle
!  Predict positions
   call pimd_coord_transform(xcart_prev,1,natom,pitransform,trotter)
   call pimd_predict_verlet(dtion,forces_pimd,mass,natom,trotter,&
&   xcart,xcart_next,xcart_prev)
!  Propagate the thermostat variables
   call pimd_nosehoover_propagate(dtion,dzeta,mass,natom,nnos,qmass,&
&   thermtemp,trotter,vel,zeta,zeta_next,zeta_prev,itimimage,pitransform)
!  Self-consistent loop
   temperature2=pimd_temperature(mass,vel)
   temp2_prev=temperature2; tol=one
   do while (tol>tolerance)
!    Recompute a (better) estimation of the velocity at time step t
     vel = (xcart_next - xcart_prev) / (two*dtion)
     dzeta=(zeta_next  - zeta_prev)  / (two*dtion)
     temperature2=pimd_temperature(mass,vel)
!    Reestimate the force
     call pimd_nosehoover_forces(dzeta,forces,forces_pimd,mass,natom,nnos,trotter,vel)
     call pimd_apply_constraint(pimd_param%constraint,constraint_output,forces_pimd,&
&     mass,natom,trotter,pimd_param%wtatcon,xcart)
!    Compute new positions
     call pimd_predict_verlet(dtion,forces_pimd,mass,natom,trotter,&
&     xcart,xcart_next,xcart_prev)
!    Propagate the thermostat variables
     call pimd_nosehoover_propagate(dtion,dzeta,mass,natom,nnos,qmass,&
&     thermtemp,trotter,vel,zeta,zeta_next,zeta_prev,itimimage,pitransform)
!    Compute variation of temperature (to check convergence of SC loop)
     tol=dabs(temperature2-temp2_prev)/dabs(temp2_prev)
     temp2_prev=temperature2
   end do ! End self-consistent loop

 end if ! itimimage==1

!Come back to primitive coordinates and velocities
 call pimd_coord_transform(xcart_next,-1,natom,pitransform,trotter)
 call pimd_coord_transform(xcart     ,-1,natom,pitransform,trotter)
 call pimd_coord_transform(xcart_prev,-1,natom,pitransform,trotter)

!Compute contributions to energy
 call pimd_energies(eharm,eharm2,epot,etotal,forces_orig,natom,spring_prim,trotter,xcart)

!Compute stress tensor at t from virial theorem
 call pimd_stresses(mass,natom,quantummass,stress_pimd,stressin,thermtemp,thermtemp,trotter,vel,volume,xcart)
 stress_pimd=-stress_pimd ! Translate pressure to stress

!############# Final operations ############################

!Print messages
 vel_cell=zero;prtstress=1;if (prtvolimg>=2) prtstress=0
 call pimd_print(pimd_param%constraint,constraint_output,&
& eharm,eharm2,epot,forces_pimd,inertmass,irestart,&
& itimimage,kt,natom,pimd_param%optcell,prtstress,prtvolimg,rprimd,&
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

!update thermostat variables
 dzeta = (three*zeta_next - four*zeta + zeta_prev)/(two * dtion)
 zeta_prev=zeta
 zeta=zeta_next
 pimd_param%zeta_prev(:,:,:,:)=zeta_prev(:,:,:,:)
 pimd_param%zeta(:,:,:,:)     =zeta(:,:,:,:)
 pimd_param%dzeta(:,:,:,:)    =dzeta(:,:,:,:)

!Free memory
 ABI_FREE(xcart)
 ABI_FREE(xcart_prev)
 ABI_FREE(xcart_next)
 ABI_FREE(forces_orig)
 ABI_FREE(forces_pimd)
 ABI_FREE(inertmass)
 ABI_FREE(quantummass)
 ABI_FREE(mass)
 ABI_FREE(spring)
 ABI_FREE(qmass)
 ABI_FREE(zeta_prev)
 ABI_FREE(zeta)
 ABI_FREE(zeta_next)
 ABI_FREE(dzeta)

end subroutine pimd_nosehoover_nvt
!!***

end module m_pimd_nosehoover
!!***

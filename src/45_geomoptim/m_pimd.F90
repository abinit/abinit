!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pimd
!! NAME
!!  m_pimd
!!
!! FUNCTION
!!  This module provides several routines and datatypes for the
!!  Path-Integral Molecular Dynamics (PIMD) implementation.
!!
!! COPYRIGHT
!! Copyright (C) 2010-2019 ABINIT group (GG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_pimd

 use defs_basis
 use defs_abitypes
 use m_abicore
 use m_dtset
 use m_errors
 use m_io_tools
 use m_random_zbq

 use m_numeric_tools,  only : uniformrandom
 use m_symtk,          only : matr3inv
 use m_geometry,       only : mkradim

 implicit none

 private

!public procedures
 public :: pimd_init
 public :: pimd_nullify
 public :: pimd_destroy
 public :: pimd_init_qtb
 public :: pimd_skip_qtb
 public :: pimd_print
 public :: pimd_is_restart
 public :: pimd_temperature
 public :: pimd_initvel
 public :: pimd_langevin_random
 public :: pimd_langevin_random_qtb
 public :: pimd_langevin_random_bar
 public :: pimd_langevin_random_init
 public :: pimd_energies
 public :: pimd_forces
 public :: pimd_langevin_forces
 public :: pimd_nosehoover_forces
 public :: pimd_stresses
 public :: pimd_diff_stress
 public :: pimd_predict_taylor
 public :: pimd_predict_verlet
 public :: pimd_nosehoover_propagate
 public :: pimd_coord_transform
 public :: pimd_force_transform
 public :: pimd_apply_constraint
 public :: pimd_mass_spring
!!***

!!****t* m_pimd/pimd_type
!! NAME
!! pimd_type
!!
!! FUNCTION
!! Datatype with the variables required to perform PIMD
!!
!! NOTES
!!
!! SOURCE

 type,public :: pimd_type
! Scalars
  integer  :: adpimd
  integer  :: constraint
  integer  :: irandom
  integer  :: nnos
  integer  :: ntypat
  integer  :: optcell
  integer  :: pitransform
  integer  :: traj_unit
  integer  :: use_qtb
  integer  :: qtb_file_unit
  real(dp) :: adpimd_gamma
  real(dp) :: vis
  real(dp) :: bmass
  real(dp) :: dtion
  real(dp) :: friction
! Arrays
  integer ,pointer  :: typat(:)      ! This pointer is associated with dtset%typat
  real(dp),pointer :: amu(:)         ! This pointer is associated with dtset%%amu_orig(:,1)
  real(dp),pointer :: mdtemp(:)      ! This pointer is associated with dtset%mdtemp
  real(dp),pointer :: pimass(:)      ! This pointer is associated with dtset%pimass
  real(dp),pointer :: qmass(:)       ! This pointer is associated with dtset%qmass
  real(dp),pointer :: strtarget(:)   ! This pointer is associated with dtset%strtarget
  real(dp),pointer :: wtatcon(:,:,:) ! This pointer is associated with dtset%wtatcon
  real(dp),allocatable :: zeta_prev(:,:,:,:)
  real(dp),allocatable :: zeta     (:,:,:,:)
  real(dp),allocatable :: zeta_next(:,:,:,:)
  real(dp),allocatable :: dzeta    (:,:,:,:)
 end type pimd_type
!!***

CONTAINS !===========================================================
!!***

!!****f* m_pimd/pimd_init
!! NAME
!!  pimd_init
!!
!! FUNCTION
!!  Initialize a datastructure of type pimd_type.
!!  Open file(s) related to this datastructure.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in current dataset
!!  is_master= TRUE if I am the master process (proc 0)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  pimd_param=datastructure of type pimd_type.
!!             several parameters for Path-Integral MD.
!!
!! PARENTS
!!      gstateimg
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_init(dtset,pimd_param,is_master)

 implicit none

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: is_master
 type(dataset_type),target,intent(in) :: dtset
 type(pimd_type),intent(inout) :: pimd_param
!Local variables-------------------------------
!scalars
 integer :: ierr
 character(len=200) :: msg

!************************************************************************

 call pimd_nullify(pimd_param)

 if((dtset%imgmov==9).or.(dtset%imgmov==10).or.(dtset%imgmov==13))then
   pimd_param%adpimd      = dtset%adpimd
   pimd_param%constraint  = dtset%pimd_constraint
   pimd_param%irandom     = dtset%irandom
   pimd_param%nnos        = dtset%nnos
   pimd_param%ntypat      = dtset%ntypat
   pimd_param%optcell     = dtset%optcell
   pimd_param%pitransform = dtset%pitransform
   pimd_param%adpimd_gamma= dtset%adpimd_gamma
   pimd_param%vis         = dtset%vis
   pimd_param%bmass       = dtset%bmass
   pimd_param%dtion       = dtset%dtion
   pimd_param%friction    = dtset%friction
   pimd_param%mdtemp      =>dtset%mdtemp
   pimd_param%pimass      =>dtset%pimass
   pimd_param%strtarget   =>dtset%strtarget
   pimd_param%amu         =>dtset%amu_orig(:,1)
   pimd_param%qmass       =>dtset%qmass
   pimd_param%typat       =>dtset%typat
   pimd_param%wtatcon     =>dtset%wtatcon
   if(dtset%imgmov==10)then
     pimd_param%use_qtb=1
     if(is_master)then
       call pimd_init_qtb(dtset,pimd_param%qtb_file_unit)
     end if
   end if
   if(dtset%imgmov==13)then
     ABI_ALLOCATE(pimd_param%zeta_prev,(3,dtset%natom,dtset%nimage,dtset%nnos))
     ABI_ALLOCATE(pimd_param%zeta     ,(3,dtset%natom,dtset%nimage,dtset%nnos))
     ABI_ALLOCATE(pimd_param%zeta_next,(3,dtset%natom,dtset%nimage,dtset%nnos))
     ABI_ALLOCATE(pimd_param%dzeta    ,(3,dtset%natom,dtset%nimage,dtset%nnos))
     pimd_param%zeta_prev=zero
     pimd_param%zeta     =zero
     pimd_param%zeta_next=zero
     pimd_param%dzeta    =zero
   end if
   if(dtset%useria==37)then
     ierr=open_file('pimd_traj.dat',msg,newunit=pimd_param%traj_unit,form='unformatted')
     if (ierr/=0) then
       MSG_ERROR(msg)
     end if
   end if
 end if

end subroutine pimd_init
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_nullify
!! NAME
!!  pimd_nullify
!!
!! FUNCTION
!!  Nullify the content of a datastructure of type pimd_type.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  pimd_param=datastructure of type pimd_type.
!!             several parameters for Path-Integral MD.
!!
!! PARENTS
!!      m_pimd
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_nullify(pimd_param)

 implicit none

!Arguments ------------------------------------
!scalars
 type(pimd_type),intent(inout) :: pimd_param

!************************************************************************

 pimd_param%adpimd       =  0
 pimd_param%constraint   =  0
 pimd_param%irandom      = -1
 pimd_param%nnos         = -1
 pimd_param%ntypat       = -1
 pimd_param%optcell      = -1
 pimd_param%pitransform  = -1
 pimd_param%qtb_file_unit= -1
 pimd_param%traj_unit    = -1
 pimd_param%use_qtb      =  0
 pimd_param%adpimd_gamma = one
 pimd_param%vis          = zero
 pimd_param%bmass        = zero
 pimd_param%dtion        = zero
 pimd_param%friction     = zero
 nullify(pimd_param%mdtemp)
 nullify(pimd_param%pimass)
 nullify(pimd_param%strtarget)
 nullify(pimd_param%amu)
 nullify(pimd_param%qmass)
 nullify(pimd_param%typat)
 nullify(pimd_param%wtatcon)

end subroutine pimd_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_destroy
!! NAME
!!  pimd_destroy
!!
!! FUNCTION
!!  Destroy the content of a datastructure of type pimd_type.
!!  Close open file(s) related to this datastructure.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  pimd_param=datastructure of type pimd_type.
!!            several parameters for PIMD.
!!
!! PARENTS
!!      gstateimg
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_destroy(pimd_param)

 implicit none

!Arguments ------------------------------------
!scalars
 type(pimd_type),intent(inout) :: pimd_param
!arrays
!Local variables-------------------------------
!scalars
 integer :: ierr
 character(len=100) :: msg
!arrays

!************************************************************************

 ABI_SFREE(pimd_param%zeta_prev)
 ABI_SFREE(pimd_param%zeta)
 ABI_SFREE(pimd_param%zeta_next)
 ABI_SFREE(pimd_param%dzeta)

 if (pimd_param%qtb_file_unit>0) then
   if (is_open(pimd_param%qtb_file_unit)) then
     ierr=close_unit(pimd_param%qtb_file_unit,msg)
   end if
 end if

 if (pimd_param%traj_unit>0) then
   if (is_open(pimd_param%traj_unit)) then
     ierr=close_unit(pimd_param%traj_unit,msg)
   end if
 end if

 call pimd_nullify(pimd_param)

end subroutine pimd_destroy
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_init_qtb
!! NAME
!!  pimd_init_qtb
!!
!! FUNCTION
!!  Only relevant for PIMD + Quantum Thermal Bath (QTB);
!!  Initialize reading of PIQTB random force file.
!!  This routine should be called only by master proc.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in current dataset
!!
!! OUTPUT
!!  qtb_file_unit=if a PIQTB_force file exists, return its file unit.
!!
!! PARENTS
!!      m_pimd
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_init_qtb(dtset,qtb_file_unit)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: qtb_file_unit
 type(dataset_type),target,intent(in) :: dtset
!Local variables-------------------------------
!scalars
 integer :: ierr,ndof_qtb,ntimimage_qtb,nimage_qtb
 real(dp) :: dtion_qtb,mdtemp_qtb
 character(len=200) :: msg

!************************************************************************

!Try to open PIQTB random force file
 ierr=open_file('piqtb_force',msg,newunit=qtb_file_unit,&
&               form='unformatted',status='old')

!Read first line of the file
 read(qtb_file_unit) dtion_qtb,ntimimage_qtb,mdtemp_qtb,nimage_qtb,ndof_qtb

!Check consistency of the read parameters with ABINIT input file
 if (abs(dtion_qtb-dtset%dtion)>tol6) then
   msg='dtion read from piqtb_force file different from dtion in input file!'
   MSG_ERROR(msg)
 end if
 if (abs(mdtemp_qtb-dtset%mdtemp(2))>tol6) then
   msg='mdtemp read from piqtb_force file different from mdtemp(2) in input file!'
   MSG_ERROR(msg)
 end if
 if (ntimimage_qtb<dtset%ntimimage) then
   msg='ntimimage read from piqtb_force file smaller than ntimimage in input file!'
   MSG_ERROR(msg)
 end if
 if (nimage_qtb/=dtset%nimage) then
   msg='nimage read from piqtb_force file different from nimage in input file!'
   MSG_ERROR(msg)
 end if
 if (ndof_qtb/=3*dtset%natom*dtset%nimage) then
   msg='Nb of degrees of freedom read from piqtb_force not consistent with input file!'
   MSG_ERROR(msg)
 end if

end subroutine pimd_init_qtb
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_skip_qtb
!! NAME
!!  pimd_skip_qtb
!!
!! FUNCTION
!!  Only relevant in case of PI-QTB:
!!  Skip a line in a QTB random force file
!!
!! INPUTS
!!  pimd_param=datastructure of type pimd_type.
!!             several parameters for Path-Integral MD.
!!
!! OUTPUT
!!
!! PARENTS
!!      gstateimg
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_skip_qtb(pimd_param)

!Arguments ------------------------------------
!scalars
 type(pimd_type),intent(in) :: pimd_param
!arrays
!Local variables-------------------------------
!scalars
 character(len=500) :: msg
!arrays

!************************************************************************

 if (pimd_param%use_qtb==0) return

 if (pimd_param%qtb_file_unit<0) then
   msg='QTB forces file unit should be positive!'
   MSG_BUG(msg)
 end if

!Skip one line QTB random forces file
 read(pimd_param%qtb_file_unit)

end subroutine pimd_skip_qtb
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_is_restart
!! NAME
!!  pimd_is_restart
!!
!! FUNCTION
!!  Determine whether this is a PIMD restart or not:
!!  test on value of velocities and corresponding temperature
!!
!! INPUTS
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  vel(3,natom,nimage)=velocities for each image of the cell
!!
!! OUTPUT
!!  pimd_is_restart=1 if temperature is not zero
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function pimd_is_restart(mass,vel,vel_cell)

!Arguments ------------------------------------
!scalars
 integer :: pimd_is_restart
!arrays
 real(dp),intent(in) :: mass(:,:),vel(:,:,:)
 real(dp),intent(in),optional :: vel_cell(:,:,:)
!Local variables-------------------------------
!scalars
 real(dp),parameter :: zero_temp=tol7
!arrays

!************************************************************************

 pimd_is_restart=0
 if (pimd_temperature(mass,vel)>zero_temp) pimd_is_restart=1
 if (present(vel_cell)) then
   if (maxval(vel_cell)>zero_temp) pimd_is_restart=pimd_is_restart+10
 end if

end function pimd_is_restart
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_temperature
!! NAME
!!  pimd_temperature
!!
!! FUNCTION
!!  Compute temperature from velocities and masses
!!
!! INPUTS
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  vel(3,natom,nimage)=velocities for each image of the cell
!!
!! OUTPUT
!!  pimd_temperature=temperature (from all images of the cell)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function pimd_temperature(mass,vel)

!Arguments ------------------------------------
!scalars
 real(dp) :: pimd_temperature
!arrays
 real(dp),intent(in) :: mass(:,:),vel(:,:,:)
!Local variables-------------------------------
!scalars
 integer :: iatom,idir,iimage,imass,natom,natom_mass,ndir,nimage,nmass
 real(dp) :: v2
 character(len=500) :: msg
!arrays

!************************************************************************

 ndir=size(vel,1);natom=size(vel,2);nimage=size(vel,3)
 natom_mass=size(mass,1);nmass=size(mass,2)
 if (ndir/=3.or.natom<=0.or.nimage<=0) then
   msg='Wrong sizes for vel array !'
   MSG_BUG(msg)
 end if
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=nimage)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

 v2=zero
 do iimage=1,nimage
   imass=min(nmass,iimage)
   do iatom=1,natom
     do idir=1,3
       v2=v2+vel(idir,iatom,iimage)*vel(idir,iatom,iimage)*mass(iatom,imass)
     end do
   end do
 end do
 pimd_temperature=v2/(dble(3*natom*nimage)*kb_HaK)

end function pimd_temperature
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_print
!! NAME
!!  pimd_print
!!
!! FUNCTION
!!  Print out results related to PIMD (for given time step)
!!
!! INPUTS
!!  constraint=type of constraint eventually applied (on a reaction coordinate)
!!  constraint_output(2)=several (real) data to be output when a constraint has been applied
!!  eharm=harmonic energy
!!  eharm_virial=harmonic energy from virial
!!  epot=potential energy
!!  forces(3,natom,trotter)=forces on atoms in each cell
!!  inertmass(natom)=inertial masses of atoms
!!  irestart=1 if this is a restart
!!  itimimage=index of time step
!!  kt=Thermal energy K_b.T
!!  natom=number of atoms
!!  optcell=option for cell evolution
!!  prtstress=flag for stress tensor printing
!!  prtvolimg=printing volume
!!  rprim(3,3)=dimensionless real space primitive translations
!!  stress(3,3,3)=stress tensor (first dimension corresponds to 3 different estimators)
!!  temperature1,temperature2=temperatures at t and t+dt
!!  traj_unit=flag activating printing of the trajectory in an external file\
!!            if >0, indicates the unit number of the trajectory file
!!  trotter=Trotter number
!!  vel(3,natom,trotter)=velocities of atoms in in each cell
!!  vel_cell(3,3)= velocities of cell parameters (time derivative of rprimd)
!!  xcart(3,natom,trotter)=cartesian coordinates of atoms in in each cell
!!  xred(3,natom,trotter)=reduced coordinates of atoms in in each cell
!!
!! OUTPUT
!!  -- only printing --
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt,pimd_nosehoover_npt
!!      pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_print(constraint,constraint_output,eharm,eharm_virial,epot,&
&          forces,inertmass,irestart,itimimage,kt,natom,optcell,prtstress,&
&          prtvolimg,rprimd,stress,temperature1,temperature2,traj_unit,&
&          trotter,vel,vel_cell,xcart,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: constraint,irestart,itimimage,natom,optcell
 integer,intent(in) :: prtstress,prtvolimg,traj_unit,trotter
 real(dp),intent(in) :: eharm,eharm_virial,epot,kt,temperature1,temperature2
!arrays
 real(dp),intent(in) :: constraint_output(2)
 real(dp),intent(in) :: forces(3,natom,trotter),inertmass(natom)
 real(dp),intent(in) :: rprimd(3,3),stress(3,3,3),vel(3,natom,trotter),vel_cell(3,3)
 real(dp),intent(in) :: xcart(3,natom,trotter),xred(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage
 real(dp) :: mtot
 character(len=500) :: msg
!arrays
 real(dp) :: acell(3),cdm(3),forcetot(3),rprim(3,3)
 real(dp),allocatable :: centroid(:,:),qudeloc(:)

!************************************************************************

!Temperature
 if(itimimage==1)then
   if(irestart==0) then
     write(msg,'(2a)') ch10,' This is a PIMD calculation from scratch'
   else if (irestart==1) then
     write(msg,'(2a)') ch10,' This is a RESTART calculation'
   end if
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a,f12.5,a)') &
&   ' In the initial configuration, the temperature is ',temperature1,' K'
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end if
 write(msg,'(2a,i5,a,f12.5,a)') ch10,&
&  ' At PIMD time step ',itimimage,', the temperature is',temperature2,' K'
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

!Energies
 write(msg,'(4a,f18.9,a,a,a,f18.9,a,a)') ch10,&
&  ' Energy:',ch10, &
&  '   Internal energy (PRIMITIVE estimator) =',onehalf*dble(natom*trotter)*kt-eharm+epot ,' Ha',ch10, &
&  '   Internal energy (VIRIAL    estimator) =',onehalf*dble(natom)*kt+eharm_virial+epot,' Ha',ch10
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

!Stress tensor and pressure
 write(msg,'(2a,3(2a,3f18.9))') ch10,&
&   ' Stress tensor from PRIMITIVE estimator (Ha/Bohr^3):',ch10, &
&   '   ',stress(1,1,1),stress(1,1,2),stress(1,1,3),ch10, &
&   '   ',stress(1,2,1),stress(1,2,2),stress(1,2,3),ch10, &
&   '   ',stress(1,3,1),stress(1,3,2),stress(1,3,3)
 if (prtstress==1) then
   call wrtout(ab_out,msg,'COLL')
 end if
 call wrtout(std_out,msg,'COLL')
 write(msg,'(a,f18.9,a)') ' Pressure (primitive estimator) =', &
&  -third*(stress(1,1,1)+stress(1,2,2)+stress(1,3,3))*HaBohr3_GPa,' GPa'
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

!Data related to constraint eventually applied
 if (constraint/=0) then
   if (constraint==1) write(msg,'(2a)') ch10,' Blue Moon Ensemble method is activated:'
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a,f18.10,2a,f18.10)') &
&     '  - Reaction coordinate =',constraint_output(1),ch10,&
&     '  - Instantaneous force on the reaction coord. =',constraint_output(2)
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end if

!Total force
 if (prtvolimg<=1) then
   forcetot=zero
   do iimage=1,trotter
     do iatom=1,natom
       do ii=1,3
         forcetot(ii)=forcetot(ii)+forces(ii,iatom,iimage)
       end do
     end do
   end do
   write(msg,'(2a,3f18.10)') ch10,' Total force=',forcetot(1:3)
   call wrtout(std_out,msg,'COLL')
 end if

!position of mass center
 if (prtvolimg<=1) then
   mtot=zero;cdm=zero
   do iimage=1,trotter
     do iatom=1,natom
       cdm(:)=cdm(:)+inertmass(iatom)*xcart(:,iatom,iimage)
       mtot=mtot+inertmass(iatom)
     end do
   end do
   cdm=cdm/mtot
   write(msg,'(3a,3f18.10)') ch10,' Center of mass:',ch10,cdm(:)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if

!Positions
 write(msg,'(2a)') ch10,' Atomic positions:'
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 do iimage=1,trotter
   select case(iimage)
   case(1)
    write(msg,'(a)') ' xred'
   case(2,3,4,5,6,7,8,9)
     write(msg,'(a,i1,a)') ' xred_',iimage,'img'
   case(10:99)
     write(msg,'(a,i2,a)') ' xred_',iimage,'img'
   case default
     write(msg,'(a,i3,a)') ' xred_',iimage,'img'
   end select
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   if (traj_unit>0) then
     call wrtout(traj_unit,msg,'COLL')
   end if
   do iatom=1,natom
     write(msg,'(3f18.10)') xred(1:3,iatom,iimage)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
     if (traj_unit>0) then
       call wrtout(traj_unit,msg,'COLL')
     end if
   end do
 end do

!Velocities
 write(msg,'(2a)') ch10,' Velocities:'
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 do iimage=1,trotter
   select case(iimage)
   case(1)
     write(msg,'(a)') ' vel'
   case(2,3,4,5,6,7,8,9)
     write(msg,'(a,i1,a)') ' vel_',iimage,'img'
   case(10:99)
     write(msg,'(a,i2,a)') ' vel_',iimage,'img'
   case default
     write(msg,'(a,i3,a)') ' vel_',iimage,'img'
   end select
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   if (traj_unit>0) then
     call wrtout(traj_unit,msg,'COLL')
   end if
   do iatom=1,natom
     write(msg,'(3f18.10)') vel(1:3,iatom,iimage)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
     if (traj_unit>0) then
       call wrtout(traj_unit,msg,'COLL')
     end if
   end do
 end do

 if (optcell>0) then

   call mkradim(acell,rprim,rprimd)

!  Time derivative of rprimd
   write(msg,'(2a)') ch10,' Time derivative of rprimd:'
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(2a,3(3f18.10))') ' vel_cell',ch10,vel_cell(:,:)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   if (traj_unit>0) then
     call wrtout(traj_unit,msg,'COLL')
   end if

!  rprimd
   write(msg,'(2a)') ch10,' Cell parameters:'
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(2a,3(3f18.10),3a,3f18.10)') ' rprim',ch10,rprim(:,:),ch10,&
&                                          ' acell',ch10,acell(:)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   if (traj_unit>0) then
     call wrtout(traj_unit,msg,'COLL')
   end if

 end if

!Centroids and wave-packet spatial spreads
 ABI_ALLOCATE(centroid,(3,natom))
 ABI_ALLOCATE(qudeloc,(natom))
 centroid=zero;qudeloc=zero
 do iimage=1,trotter
   do iatom=1,natom
     do ii=1,3
       centroid(ii,iatom)=centroid(ii,iatom)+xcart(ii,iatom,iimage)
     end do
   end do
 end do
 centroid=centroid/dble(trotter)
 do iimage=1,trotter
   do iatom=1,natom
     do ii=1,3
       qudeloc(iatom)=qudeloc(iatom)+((xcart(ii,iatom,iimage)-centroid(ii,iatom))**2)
     end do
   end do
 end do
 qudeloc(:)=sqrt(qudeloc(:)/dble(trotter))
 write(msg,'(4a)') ch10,' Centroids and wave-packet spatial spreads (cart. coord.):',ch10,&
&  ' iat        centroid_x        centroid_y        centroid_z    spatial_spread'
 call wrtout(std_out,msg,'COLL')
 do iatom=1,natom
   write(msg,'(i4,4f18.10)') iatom,centroid(1:3,iatom),qudeloc(iatom)
   call wrtout(std_out,msg,'COLL')
 end do
 ABI_DEALLOCATE(centroid)
 ABI_DEALLOCATE(qudeloc)

!Fake statement
 return;ii=prtvolimg

end subroutine pimd_print
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_initvel
!! NAME
!!  pimd_initvel
!!
!! FUNCTION
!!  Initialize velocities for PIMD with a gaussian distribution
!!  fixing the center of mass
!!  and eventually applying a constraint on atomic positions
!!
!! INPUTS
!!  constraint=type of constraint to be applied
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  temperature=temperature used to define velocities
!!  trotter=Trotter number
!!  wtatcon(3,natom)=weights for atomic constraints
!!
!! OUTPUT
!!  vel(3,natom,trotter)=velocities of atoms in in each cell
!!
!! SIDE EFFECTS
!!  iseed=seed for random number generator
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt,pimd_nosehoover_npt
!!      pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_initvel(iseed,mass,natom,temperature,trotter,vel,constraint,wtatcon)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: constraint,natom,trotter
 integer,intent(inout) :: iseed
 real(dp),intent(in) :: temperature
!arrays
 real(dp),intent(in) :: mass(:,:),wtatcon(3,natom)
 real(dp),intent(out) :: vel(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,imass,natom_mass,nmass
 real(dp) :: mtot,rescale_vel
 character(len=500) :: msg
!arrays
 real(dp) :: mvini(3)

!************************************************************************

 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

!Compute total mass (of non-constrained atoms)
 if (constraint==0) then
   mtot=sum(mass(1:natom,1:nmass))
 else
   mtot=zero
   do iatom=1,natom
     if (all(abs(wtatcon(:,iatom))<tol8)) mtot=mtot+sum(mass(iatom,1:nmass))
   end do
 end if
 if (nmass==1) mtot=mtot*dble(trotter)

!Initialize randomly the velocities
 do iimage=1,trotter
   imass=min(nmass,iimage)
   do iatom=1,natom
     do ii=1,3
       vel(ii,iatom,iimage)=sqrt(kb_HaK*temperature/mass(iatom,imass))*cos(two_pi*uniformrandom(iseed))
       vel(ii,iatom,iimage)=vel(ii,iatom,iimage)*sqrt(-two*log(uniformrandom(iseed)))
     end do
   end do
 end do

!Cancel velocities of constrained atoms
 if (constraint/=0) then
   do iimage=1,trotter
     do iatom=1,natom
       if (any(abs(wtatcon(:,iatom))>=tol8)) vel(:,iatom,iimage)=zero
     end do
   end do
 end if

!Make sure that the (sum of m_i v_i) at step zero is zero
 mvini=zero
 do iimage=1,trotter
   imass=min(nmass,iimage)
   do iatom=1,natom
     do ii=1,3
      mvini(ii)=mvini(ii)+mass(iatom,imass)*vel(ii,iatom,iimage)
     end do
   end do
 end do
 if (constraint==0) then
   do iimage=1,trotter
     do iatom=1,natom
       do ii=1,3
         vel(ii,iatom,iimage)=vel(ii,iatom,iimage)-(mvini(ii)/mtot)
       end do
     end do
   end do
 else
   do iimage=1,trotter
     do iatom=1,natom
       if (all(abs(wtatcon(:,iatom))<tol8)) vel(:,iatom,iimage)=vel(:,iatom,iimage)-(mvini(:)/mtot)
     end do
   end do
 end if

!Now rescale the velocities to give the exact temperature
 rescale_vel=sqrt(temperature/pimd_temperature(mass,vel))
 vel(:,:,:)=vel(:,:,:)*rescale_vel

end subroutine pimd_initvel
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_langevin_random
!! NAME
!!  pimd_langevin_random
!!
!! FUNCTION
!!  Generate a set of random numbers to be used for PIMD Langevin algorithm
!!
!! INPUTS
!!  irandom=option for random number generator:
!!          1:uniform random routine provided within Abinit package
!!          2:Fortran 90 random number generator
!!          3:ZBQLU01 non deterministic random number generator
!!  langev(natom,mass_dim)=Langevin factors (mass_dim=1 or trotter)
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  trotter=Trotter number
!!  zeroforce=flag; if 1 keep sum of forces equal to zero
!!
!! OUTPUT
!!  alea(3,natom,trotter)=set of random numbers
!!
!! SIDE EFFECTS
!!  iseed=seed for random number generator (used only if irandom=1)
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_langevin_random(alea,irandom,iseed,langev,mass,natom,trotter,zeroforce)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: irandom,natom,trotter,zeroforce
 integer,intent(inout) :: iseed
!arrays
 real(dp),intent(in) :: langev(:,:),mass(:,:)
 real(dp),intent(out) :: alea(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,imass,nmass,natom_mass
 real(dp) :: mtot,r1,r2
 character(len=500) :: msg
!arrays
 real(dp) :: total(3)

!************************************************************************

 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

 mtot=sum(mass(1:natom,1:nmass))
 if (nmass==1) mtot=mtot*trotter

!Draw random numbers
 do iimage=1,trotter
   do iatom=1,natom
     do ii=1,3
       select case(irandom)
       case(1)
         r1=uniformrandom(iseed)
         r2=uniformrandom(iseed)
       case(2)
         call random_number(r1)
         call random_number(r2)
       case(3)
         r1=ZBQLU01(zero)
         r2=ZBQLU01(zero)
       end select
       alea(ii,iatom,iimage)= cos(two_pi*r1)*sqrt(-log(r2)*two)
     end do
   end do
 end do

!Make sure that the sum of random forces is zero
 if(zeroforce==1)then
   total=zero
   mtot=zero
   do iimage=1,trotter
     imass=min(nmass,iimage)
     do iatom=1,natom
       mtot=mtot+mass(iatom,imass)
     end do
   end do
   do iimage=1,trotter
     imass=min(nmass,iimage)
     do iatom=1,natom
       do ii=1,3
         total(ii)=total(ii)+langev(iatom,imass)*alea(ii,iatom,iimage)
       end do
     end do
   end do
   do iimage=1,trotter
     imass=min(nmass,iimage)
     do iatom=1,natom
       do ii=1,3
         alea(ii,iatom,iimage)= alea(ii,iatom,iimage)- &
&        (total(ii)*mass(iatom,imass))/(langev(iatom,imass)*mtot)
       end do
     end do
   end do
!  now random forces have been rescaled so that their sum is zero
 end if

end subroutine pimd_langevin_random
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_langevin_random_qtb
!! NAME
!!  pimd_langevin_random_qtb
!!
!! FUNCTION
!!  Read a set of random forces (atm units) to be used for PIMD QTB algorithm
!!
!! INPUTS
!!  langev(natom,mass_dim)=Langevin factors (mass_dim=1 or trotter)
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  qtb_file_unit= ramdom forces file unit
!!  trotter=Trotter number
!!  zeroforce=flag; if 1 keep sum of forces equal to zero
!!
!! OUTPUT
!!  alea(3,natom,trotter)=set of random forces
!!
!! PARENTS
!!      pimd_langevin_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_langevin_random_qtb(alea,langev,mass,natom,qtb_file_unit,trotter,zeroforce)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,qtb_file_unit,trotter,zeroforce
!arrays
 real(dp),intent(in) :: langev(:,:),mass(:,:)
 real(dp),intent(out) :: alea(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,imass,nmass,natom_mass
 real(dp) :: mtot
 character(len=500) :: msg
!arrays
 real(sp) :: alea_sp(3,natom,trotter)
 real(dp) :: total(3)

!************************************************************************

 if (qtb_file_unit<0) then
   msg='QTB forces file unit should be positive!'
   MSG_BUG(msg)
 end if
 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

!Read QTB random forces
 read(qtb_file_unit) alea_sp(1:3,1:natom,1:trotter)
 alea(:,:,:)=dble(alea_sp(:,:,:))

!Make sure that the sum of random forces is zero
 if(zeroforce==1)then
   total=zero
   mtot=zero
   do iimage=1,trotter
     imass=min(nmass,iimage)
     do iatom=1,natom
       mtot=mtot+mass(iatom,imass)
     end do
   end do
   do iimage=1,trotter
     imass=min(nmass,iimage)
     do iatom=1,natom
       do ii=1,3
         total(ii)=total(ii)+langev(iatom,imass)*alea(ii,iatom,iimage)
       end do
     end do
   end do
   do iimage=1,trotter
     imass=min(nmass,iimage)
     do iatom=1,natom
       do ii=1,3
         alea(ii,iatom,iimage)= alea(ii,iatom,iimage)- &
&        (total(ii)*mass(iatom,imass))/(langev(iatom,imass)*mtot)
       end do
     end do
   end do
!  now random forces have been rescaled so that their sum is zero
 end if

end subroutine pimd_langevin_random_qtb
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_langevin_random_bar
!! NAME
!!  pimd_langevin_random_bar
!!
!! FUNCTION
!!  Generate a set of random numbers to be used for the barostat of PIMD Langevin algorithm
!!
!! INPUTS
!!  irandom=option for random number generator:
!!          1:uniform random routine provided within Abinit package
!!          2:Fortran 90 random number generator
!!          3:ZBQLU01 non deterministic random number generator
!!
!! OUTPUT
!!  alea_bar(3,3)=set of random numbers
!!
!! SIDE EFFECTS
!!  iseed=seed for random number generator (used only if irandom=1)
!!
!! PARENTS
!!      pimd_langevin_npt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_langevin_random_bar(alea_bar,irandom,iseed)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: irandom
 integer,intent(inout) :: iseed
!arrays
 real(dp),intent(out) :: alea_bar(3,3)
!Local variables-------------------------------
!scalars
 integer :: ii,jj
 real(dp) :: r1,r2
!arrays

!************************************************************************

!Draw random numbers
 do ii=1,3
   do jj=1,3
     select case(irandom)
     case(1)
       r1=uniformrandom(iseed)
       r2=uniformrandom(iseed)
     case(2)
       call random_number(r1)
       call random_number(r2)
     case(3)
       r1=ZBQLU01(zero)
       r2=ZBQLU01(zero)
     end select
     alea_bar(ii,jj)= cos(two_pi*r1)*sqrt(-log(r2)*two)
   end do
 end do

!Symmetrize
 alea_bar(1,2)=half*(alea_bar(1,2)+alea_bar(2,1))
 alea_bar(1,3)=half*(alea_bar(1,3)+alea_bar(3,1))
 alea_bar(2,3)=half*(alea_bar(2,3)+alea_bar(3,2))
 alea_bar(2,1)=alea_bar(1,2)
 alea_bar(3,1)=alea_bar(1,3)
 alea_bar(3,2)=alea_bar(2,3)

end subroutine pimd_langevin_random_bar
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_langevin_random_init
!! NAME
!!  pimd_langevin_random_init
!!
!! FUNCTION
!!  Initialize random number generator to be used for PIMD Langevin algorithm
!!
!! INPUTS
!!  irandom=option for random number generator:
!!          1:uniform random routine provided within Abinit package
!!          2:Fortran 90 random number generator
!!          3:ZBQLU01 non deterministic random number generator
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  iseed=seed for random number generator (used only if irandom=1)
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_langevin_random_init(irandom,iseed)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: irandom
 integer,intent(inout) :: iseed

!************************************************************************

 if (irandom==3) then
   call ZBQLINI(0)
 end if

!Fake statement
 return;if (.false.) iseed=zero

end subroutine pimd_langevin_random_init
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_energies
!! NAME
!!  pimd_energies
!!
!! FUNCTION
!!  In the case od PIMD, compute the several contribution to total energy
!!
!! INPUTS
!!  etotal_img(trotter)= energy (from DFT) for each cell
!!  forces(3,natom,trotter)=forces (from DFT) on atoms in each cell
!!  natom=number of atoms
!!  spring(natom)=spring constants in the primitive scheme
!!  trotter=Trotter number
!!  xcart(3,natom,trotter)=cartesian coordinates of atoms in each cell at t
!!
!! OUTPUT
!!  eharm       =harmonic energy
!!  eharm_virial=harmonic energy from virial estimator
!!  epot        =potential energy
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt,pimd_nosehoover_npt
!!      pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_energies(eharm,eharm_virial,epot,etotal_img,forces,natom,spring,trotter,xcart)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,trotter
 real(dp),intent(out) :: eharm,eharm_virial,epot
!arrays
 real(dp),intent(in) :: etotal_img(trotter),forces(3,natom,trotter)
 real(dp),intent(in) :: xcart(3,natom,trotter)
 real(dp),intent(in) :: spring(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,iimagep
!arrays
 real(dp),allocatable :: centroid(:,:)

!************************************************************************

!Compute the centroid
 ABI_ALLOCATE(centroid,(3,natom))
 centroid=zero
 do iimage=1,trotter
   do iatom=1,natom
     do ii=1,3
       centroid(ii,iatom)=centroid(ii,iatom)+xcart(ii,iatom,iimage)
     end do
   end do
 end do
 centroid=centroid/dble(trotter)

!Potential energy
 epot=sum(etotal_img(1:trotter))/dble(trotter)

!Harmonic energy
 eharm=zero
 do iimage=1,trotter
   iimagep=iimage+1;if(iimage==trotter)iimagep=1
   do iatom=1,natom
     do ii=1,3
       eharm=eharm+half*spring(iatom)*((xcart(ii,iatom,iimagep)-xcart(ii,iatom,iimage))**2)
     end do
   end do
 end do

!Harmonic energy from virial estimator
 eharm_virial=zero
 do iimage=1,trotter
   do iatom=1,natom
     do ii=1,3
       eharm_virial=eharm_virial-(xcart(ii,iatom,iimage)-centroid(ii,iatom)) &
&              *forces(ii,iatom,iimage)
     end do
   end do
 end do
 eharm_virial=eharm_virial/dble(two*trotter)

 ABI_DEALLOCATE(centroid)

end subroutine pimd_energies
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_forces
!! NAME
!!  pimd_forces
!!
!! FUNCTION
!!  Modify forces in order to take into account PIMD contribution
!!
!! INPUTS
!!  natom=number of atoms
!!  spring(natom,spring_dim)=spring constants (spring_dim=1 or trotter)
!!  transform=coordinate transformation:
!!            0: no tranformation
!!            1: normal mode transformation
!!            2: staging transformation
!!  trotter=Trotter number
!!  xcart(3,natom,trotter)=cartesian coordinates of atoms in each cell
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  forces(3,natom,trotter)=
!!    at input:  forces from electronic calculation
!!    at output: forces from electronic calculation + quantum spring contribution
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt,pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_forces(forces,natom,spring,transform,trotter,xcart)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,transform,trotter
!arrays
 real(dp),intent(in) :: xcart(3,natom,trotter)
 real(dp),intent(in) :: spring(:,:)
 real(dp),intent(inout) :: forces(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,iimagem,iimagep,ispring,natom_spring,nspring
 character(len=500) :: msg
!arrays

!************************************************************************

 natom_spring=size(spring,1);nspring=size(spring,2)
 if (natom/=natom_spring.or.(nspring/=1.and.nspring/=trotter)) then
   msg='Wrong dimensions for array spring !'
   MSG_BUG(msg)
 end if

 if (transform==0) then
   do iimage=1,trotter
     ispring=min(nspring,iimage)
     iimagep=iimage+1; iimagem=iimage-1
     if(iimage==trotter) iimagep=1
     if(iimage==1)       iimagem=trotter
     do iatom=1,natom
       do ii=1,3
         forces(ii,iatom,iimage)= &
&               forces(ii,iatom,iimage)/dble(trotter) &
&             - spring(iatom,ispring)*(two*xcart(ii,iatom,iimage)-xcart(ii,iatom,iimagem) &
&                                                                -xcart(ii,iatom,iimagep))
       end do
     end do
   end do

 else
   do iimage=1,trotter
     ispring=min(nspring,iimage)
     do iatom=1,natom
       do ii=1,3
         forces(ii,iatom,iimage)= &
&               forces(ii,iatom,iimage)/dble(trotter) &
&             - spring(iatom,ispring)*xcart(ii,iatom,iimage)
       end do
     end do
   end do

 end if

end subroutine pimd_forces
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_langevin_forces
!! NAME
!!  pimd_langevin_forces
!!
!! FUNCTION
!!  Compute Langevin contribution to PIMD forces
!!
!! INPUTS
!!  alea(3,natom,trotter)=set of random numbers
!!  forces(3,natom,trotter)=forces without Langevin contribution
!!  friction=friction factor
!!  langev(natom,mass_dim)=Langevin factors (mass_dim=1 or trotter)
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  trotter=Trotter number
!!  vel(3,natom,trotter)=velocities of atoms in each cell
!!
!! OUTPUT
!!  forces_langevin(3,natom,trotter)=forces including Langevin contribution
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_langevin_forces(alea,forces,forces_langevin,friction,&
&                               langev,mass,natom,trotter,vel)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,trotter
 real(dp),intent(in) :: friction
!arrays
 real(dp),intent(in) :: alea(3,natom,trotter),forces(3,natom,trotter)
 real(dp),intent(in) :: langev(:,:),mass(:,:),vel(3,natom,trotter)
 real(dp),intent(out) :: forces_langevin(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,imass,natom_mass,nmass
 character(len=500) :: msg
!arrays

!************************************************************************

 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

 do iimage=1,trotter
   imass=min(nmass,iimage)
   do iatom=1,natom
     do ii=1,3
       forces_langevin(ii,iatom,iimage)=forces(ii,iatom,iimage) &
&                    + langev(iatom,imass)*alea(ii,iatom,iimage) &
&                    - friction*mass(iatom,imass)*vel(ii,iatom,iimage)
     end do
   end do
 end do

end subroutine pimd_langevin_forces
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_noseehoover_forces
!! NAME
!!  pimd_nosehoover_forces
!!
!! FUNCTION
!!  Compute Nose-Hoover contribution to PIMD forces
!!  by adding friction force of thermostat number one
!!
!! INPUTS
!!  dzeta(3,natom,trotter,nnos)=variables of thermostats, in (atomic time unit)^(-1)
!!                              used only when a coordinate transformation is applied (transfom/=0)
!!  forces(3,natom,trotter)=forces without Nose-Hoover contribution
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  nnos=number of thermostats
!!  trotter=Trotter number
!!  vel(3,natom,trotter)=velocities of atoms in each cell
!!
!! OUTPUT
!!  forces_nosehoover(3,natom,trotter)=forces including thermostat contribution
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_nosehoover_forces(dzeta,forces,forces_nosehoover,mass,natom,&
&                                 nnos,trotter,vel)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nnos,trotter
!arrays
 real(dp),intent(in) :: dzeta(3,natom,trotter,nnos),forces(3,natom,trotter)
 real(dp),intent(in) :: vel(3,natom,trotter)
 real(dp),intent(in) :: mass(:,:)
 real(dp),intent(out) :: forces_nosehoover(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,imass,natom_mass,nmass
 character(len=500) :: msg
!arrays

!************************************************************************

 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

  do iimage=1,trotter
    imass=min(nmass,iimage)
    do iatom=1,natom
      do ii=1,3
        forces_nosehoover(ii,iatom,iimage)=forces(ii,iatom,iimage) &
&           - mass(iatom,imass)*dzeta(ii,iatom,iimage,1)*vel(ii,iatom,iimage)
      end do
    end do
  end do

end subroutine pimd_nosehoover_forces
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_stresses
!! NAME
!!  pimd_stresses
!!
!! FUNCTION
!!  In the case od PIMD, compute the pressure tensor from virial theorem
!!
!! INPUTS
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  quantummass(natom)=quantum masses of atoms
!!  stressin(3,3,trotter)=electronic stress tensor for each image
!!  temperature=temperature (could be instantaneous temp. or thermostat temp.)
!!  temperature_therm=thermostat temperature
!!  trotter=Trotter number
!!  vel(3,natom,trotter)=velocities of atoms in each cell
!!  volume=volume of each cell (common to all cells)
!!  xcart(3,natom,trotter)=cartesian coordinates of atoms in each cell at t
!!  temperature=thermostat temperature
!!
!! OUTPUT
!!  stress_pimd(3,3,3)=stress tensor for PIMD
!!                     First dimension (3) corresponds to 3 different pressure estimators
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt,pimd_nosehoover_npt
!!      pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_stresses(mass,natom,quantummass,stress_pimd,stressin,&
&                        temperature,temperature_therm,trotter,vel,volume,xcart)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,trotter
 real(dp),intent(in) :: temperature,temperature_therm,volume
!arrays
 real(dp),intent(in) :: vel(3,natom,trotter),xcart(3,natom,trotter)
 real(dp),intent(in) :: mass(:,:),stressin(3,3,trotter),quantummass(natom)
 real(dp),intent(out) :: stress_pimd(3,3,3)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,imass,jj,natom_mass,nmass,iimagep
 real(dp) :: stress_tmp(3,3,3),omega2,kt,kt_therm
 character(len=500) :: msg

!************************************************************************

 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

 kt=temperature*kb_HaK
 kt_therm=temperature_therm*kb_HaK
 stress_pimd=zero

!I-PRIMITIVE ESTIMATOR
!1-Kinetic part
 do ii=1,3
   stress_pimd(1,ii,ii)=dble(natom)*dble(trotter)*kt/volume
 end do

!2-Potential part
 do iimage=1,trotter
   do ii=1,3
     do jj=1,3
       stress_pimd(1,ii,jj)=stress_pimd(1,ii,jj)-stressin(ii,jj,iimage)/dble(trotter)
       !minus to convert stress into pressure
     end do
   end do
 end do

!3-Contribution from springs
 omega2=dble(trotter)*kt_therm*kt_therm
 do iimage=1,trotter
   iimagep=iimage+1
   if(iimage==trotter) then
     iimagep=1
   end if
   do iatom=1,natom
     do ii=1,3
       do jj=1,3
         stress_pimd(1,ii,jj)=stress_pimd(1,ii,jj)-quantummass(iatom)*omega2* &
&        (xcart(ii,iatom,iimagep)-xcart(ii,iatom,iimage))*  &
&        (xcart(jj,iatom,iimagep)-xcart(jj,iatom,iimage))/volume
       end do
     end do
   end do
 end do

!II-Average of classical pressures
!1-Kinetic part
 do iimage=1,trotter
  imass=min(nmass,iimage)
  do iatom=1,natom
    do ii=1,3
      do jj=1,3
        stress_pimd(2,ii,jj)=stress_pimd(2,ii,jj)+ &
&          mass(iatom,imass)*vel(ii,iatom,iimage)*vel(jj,iatom,iimage)/volume
      end do
    end do
  end do
 end do

!2-Contribution from electronic stress
 do iimage=1,trotter
   do ii=1,3
     do jj=1,3
       stress_pimd(2,ii,jj)=stress_pimd(2,ii,jj)-stressin(ii,jj,iimage)
     end do
   end do
 end do
 do ii=1,3
   do jj=1,3
     stress_pimd(2,ii,jj)=stress_pimd(2,ii,jj)/dble(trotter)
   end do
 end do

!III-pressure from VIRIAL estimator: stress_pimd(3,:,:)
!1-kinetic part
! do ii=1,3
!   stress_pimd(3,ii,ii)=dfloat(natom)*kt/volume
! end do

!Symmetrize internal pressure
 stress_tmp=stress_pimd
 stress_pimd(:,2,1)=half*(stress_tmp(:,1,2)+stress_tmp(:,2,1))
 stress_pimd(:,1,2)=half*(stress_tmp(:,1,2)+stress_tmp(:,2,1))
 stress_pimd(:,3,1)=half*(stress_tmp(:,1,3)+stress_tmp(:,3,1))
 stress_pimd(:,1,3)=half*(stress_tmp(:,1,3)+stress_tmp(:,3,1))
 stress_pimd(:,2,3)=half*(stress_tmp(:,3,2)+stress_tmp(:,2,3))
 stress_pimd(:,3,2)=half*(stress_tmp(:,3,2)+stress_tmp(:,2,3))

end subroutine pimd_stresses
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_diff_stress
!! NAME
!!  pimd_diff_stress
!!
!! FUNCTION
!!  Compute the difference between the stress tensor and the stress target
!!
!! INPUTS
!!  stress_pimd(3,3,3)=stress tensor for PIMD
!!                     Last dimension (3) corresponds to 3 different pressure estimators
!!  stress_target(6)=stress target
!!
!! OUTPUT
!!  pimd_diff_stress(3,3,3)=difference between stresses and stress target
!!                          First dimension (3) corresponds to 3 different pressure estimators
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function pimd_diff_stress(stress_pimd,stress_target)

!Arguments ------------------------------------
!scalars
!arrays
 real(dp),intent(in) :: stress_pimd(3,3,3),stress_target(6)
 real(dp) :: pimd_diff_stress(3,3)
!Local variables-------------------------------
!scalars
 integer :: ii,jj
!arrays
 real(dp) :: stress_pimd2(3,3)

!************************************************************************

!Choice: the primitive estimator for pressure is chosen
 do ii=1,3
   do jj=1,3
     stress_pimd2(ii,jj)=stress_pimd(1,ii,jj)
   end do
 end do

!+stress_target instead of - because it is translated from stress to pressure tensor
 pimd_diff_stress(1,1)=stress_pimd2(1,1)+stress_target(1)
 pimd_diff_stress(2,2)=stress_pimd2(2,2)+stress_target(2)
 pimd_diff_stress(3,3)=stress_pimd2(3,3)+stress_target(3)
 pimd_diff_stress(2,3)=stress_pimd2(2,3)+stress_target(4)
 pimd_diff_stress(3,2)=stress_pimd2(3,2)+stress_target(4)
 pimd_diff_stress(1,3)=stress_pimd2(1,3)+stress_target(5)
 pimd_diff_stress(3,1)=stress_pimd2(3,1)+stress_target(5)
 pimd_diff_stress(1,2)=stress_pimd2(1,2)+stress_target(6)
 pimd_diff_stress(2,1)=stress_pimd2(2,1)+stress_target(6)

end function pimd_diff_stress
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_predict_taylor
!! NAME
!!  pimd_predict_taylor
!!
!! FUNCTION
!!  Predict new atomic positions using a Taylor algorithm (first time step) - for PIMD
!!
!! INPUTS
!!  dtion=time step
!!  forces(3,natom,trotter)=PIMD forces on atoms in each cell
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  trotter=Trotter number
!!  vel(3,natom,trotter)=velocities of atoms in each cell
!!  xcart(3,natom,trotter)=cartesian coordinates of atoms in each cell at t
!!
!! OUTPUT
!!  xcart_next(3,natom,trotter)=cartesian coordinates of atoms in each cell at t+dt
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_langevin_nvt,pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_predict_taylor(dtion,forces,mass,natom,trotter,vel,xcart,xcart_next)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,trotter
 real(dp),intent(in) :: dtion
!arrays
 real(dp),intent(in) :: forces(3,natom,trotter),vel(3,natom,trotter),xcart(3,natom,trotter)
 real(dp),intent(in) :: mass(:,:)
 real(dp),intent(out) :: xcart_next(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,imass,natom_mass,nmass
 character(len=500) :: msg
!arrays

!************************************************************************

 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

 do iimage=1,trotter
   imass=min(nmass,iimage)
   do iatom=1,natom
     do ii=1,3
       xcart_next(ii,iatom,iimage)=xcart(ii,iatom,iimage) &
&             + half*dtion*dtion*forces(ii,iatom,iimage)/mass(iatom,imass) &
&             + dtion*vel(ii,iatom,iimage)
     end do
   end do
 end do

end subroutine pimd_predict_taylor
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_predict_verlet
!! NAME
!!  pimd_predict_verlet
!!
!! FUNCTION
!!  Predict new atomic positions using a Verlet algorithm - for PIMD
!!
!! INPUTS
!!  dtion=time step
!!  forces(3,natom,trotter)=PIMD forces on atoms in each cell
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  trotter=Trotter number
!!  xcart(3,natom,trotter)=cartesian coordinates of atoms in each cell at t
!!  xcart_prev(3,natom,trotter)=cartesian coordinates of atoms in each cell at t-dt
!!
!! OUTPUT
!!  xcart_next(3,natom,trotter)=cartesian coordinates of atoms in each cell at t+dt
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_langevin_nvt,pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_predict_verlet(dtion,forces,mass,natom,trotter,xcart,xcart_next,xcart_prev)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,trotter
 real(dp),intent(in) :: dtion
!arrays
 real(dp),intent(in) :: forces(3,natom,trotter)
 real(dp),intent(in) :: xcart(3,natom,trotter),xcart_prev(3,natom,trotter)
 real(dp),intent(in) :: mass(:,:)
 real(dp),intent(out) :: xcart_next(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,imass,natom_mass,nmass
 character(len=500) :: msg
!arrays

!************************************************************************

 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

 do iimage=1,trotter
   imass=min(nmass,iimage)
   do iatom=1,natom
     do ii=1,3
       xcart_next(ii,iatom,iimage)= &
&         two*xcart(ii,iatom,iimage) &
&       - xcart_prev(ii,iatom,iimage) &
&       + dtion*dtion*forces(ii,iatom,iimage)/mass(iatom,imass)
     end do
   end do
 end do

end subroutine pimd_predict_verlet
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_nosehoover_propagate
!! NAME
!!  pimd_nosehoover_propagate
!!
!! FUNCTION
!!  Propagate thermostat variables (Nose-Hoover algorithm) - for PIMD
!!
!! INPUTS
!!  dtion=time step
!!  dzeta(3,natom,trotter,nnos)=time derivative of zeta at time t
!!  itimimage=index of time step
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  nnos=number of thermostats
!!  qmass(nnos)=masses of thermostats
!!  temperature=temperature
!!  transform=coordinate transformation:
!!            0: no tranformation
!!            1: normal mode transformation
!!            2: staging transformation
!!  trotter=Trotter number
!!  vel(3,natom,trotter)=velocities of atoms in each cell
!!  zeta(3,natom,trotter,nnos)=variables of thermostats, in (atomic time unit)^(-1)
!!             used only when no coordinate transformation is applied (transfom==0)
!!             at time t
!!  zeta_prev(3,natom,trotter,nnos)=previous value of zeta (t-dt)
!!
!! OUTPUT
!!  zeta_next(3,natom,trotter,nnos)=next value of zeta (t+dt)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_nosehoover_propagate(dtion,dzeta,mass,natom,nnos,qmass,temperature,&
&                                    trotter,vel,zeta,zeta_next,zeta_prev,itimimage,transform)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimimage,natom,nnos,transform,trotter
 real(dp),intent(in) :: dtion,temperature
!arrays
 real(dp),intent(in) :: qmass(nnos),vel(3,natom,trotter)
 real(dp),intent(in) :: mass(:,:)
 real(dp),intent(in) :: dzeta(3,natom,trotter,nnos),zeta_prev(3,natom,trotter,nnos)
 real(dp),intent(in) :: zeta(3,natom,trotter,nnos)
 real(dp),intent(out) :: zeta_next(3,natom,trotter,nnos)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,inos,natom_mass,nmass
 real(dp) :: kt
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: thermforces(:,:,:,:)

!************************************************************************

 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if
 if (nnos<3) then
   msg='Not available for nnos<3 !'
   MSG_BUG(msg)
 end if

 kt=temperature*kb_HaK
 ABI_ALLOCATE(thermforces,(3,natom,trotter,nnos))

!forces on the thermostats
 do inos=1,nnos
   if(inos==1)then
     do iimage=1,trotter
       do iatom=1,natom
         do ii=1,3
           select case(transform)
           case(0)
             thermforces(ii,iatom,iimage,1)=mass(iatom,1)*(vel(ii,iatom,iimage)**2)-kt
           case(1)
             thermforces(ii,iatom,iimage,1)=mass(iatom,iimage)*(vel(ii,iatom,iimage)**2)-kt
           case(2)
             thermforces(ii,iatom,iimage,1)=mass(iatom,iimage)*(vel(ii,iatom,iimage)**2)-kt
           end select
         end do
       end do
     end do
   elseif(inos==(nnos-1))then
     do iimage=1,trotter
       do iatom=1,natom
         do ii=1,3
           thermforces(ii,iatom,iimage,inos)=&
&          (qmass(nnos-2)*(dzeta(ii,iatom,iimage,nnos-2)**2)-kt) + &
&          (qmass(nnos  )*(dzeta(ii,iatom,iimage,nnos)  **2)-kt)
         end do
       end do
     end do
   elseif(inos==nnos)then
     do iimage=1,trotter
       do iatom=1,natom
         do ii=1,3
           thermforces(ii,iatom,iimage,inos)=&
&          qmass(nnos-1)*(dzeta(ii,iatom,iimage,nnos-1)**2)-kt
         end do
       end do
     end do
   else
     do iimage=1,trotter
       do iatom=1,natom
         do ii=1,3
           thermforces(ii,iatom,iimage,inos)=&
&          qmass(inos-1)*(dzeta(ii,iatom,iimage,inos-1)**2)-kt
         end do
       end do
     end do
   end if
 end do

 select case(itimimage)
 case(1) !taylor

 do inos=1,nnos
   if(inos==1)then
     zeta_next(:,:,:,1)=zeta(:,:,:,1)+ dzeta(:,:,:,1)*dtion + &
&       (thermforces(:,:,:,1)-qmass(1)*dzeta(:,:,:,1)*dzeta(:,:,:,2))* &
&       dtion*dtion/(two*qmass(1))
   elseif(inos==(nnos-1))then
     zeta_next(:,:,:,inos)=zeta(:,:,:,inos)+ dzeta(:,:,:,inos)*dtion + &
&       (thermforces(:,:,:,inos)-qmass(inos)*dzeta(:,:,:,inos)*dzeta(:,:,:,nnos))* &
&       dtion*dtion/(two*qmass(inos))
   elseif(inos==nnos)then
     zeta_next(:,:,:,inos)=zeta(:,:,:,inos)+ dzeta(:,:,:,inos)*dtion + &
&       (thermforces(:,:,:,inos)-qmass(inos)*dzeta(:,:,:,nnos-1)*dzeta(:,:,:,nnos))* &
&       dtion*dtion/(two*qmass(inos))
   else
     zeta_next(:,:,:,inos)=zeta(:,:,:,inos)+ dzeta(:,:,:,inos)*dtion + &
&       (thermforces(:,:,:,inos)-qmass(inos)*dzeta(:,:,:,inos)*dzeta(:,:,:,inos+1))* &
&       dtion*dtion/(two*qmass(inos))
   end if
 end do

 case default !verlet

 do inos=1,nnos
   if(inos==1)then
     zeta_next(:,:,:,1)=two*zeta(:,:,:,1) - zeta_prev(:,:,:,1) + &
&       (thermforces(:,:,:,1)-qmass(1)*dzeta(:,:,:,1)*dzeta(:,:,:,2))* &
&       dtion*dtion/qmass(1)
   elseif(inos==(nnos-1))then
     zeta_next(:,:,:,inos)=two*zeta(:,:,:,inos) - zeta_prev(:,:,:,inos) + &
&       (thermforces(:,:,:,inos)-qmass(inos)*dzeta(:,:,:,inos)*dzeta(:,:,:,nnos))* &
&       dtion*dtion/qmass(inos)
   elseif(inos==nnos)then
     zeta_next(:,:,:,inos)=two*zeta(:,:,:,inos) - zeta_prev(:,:,:,inos) + &
&       (thermforces(:,:,:,inos)-qmass(inos)*dzeta(:,:,:,nnos-1)*dzeta(:,:,:,nnos))* &
&       dtion*dtion/qmass(inos)
   else
     zeta_next(:,:,:,inos)=two*zeta(:,:,:,inos) - zeta_prev(:,:,:,inos) + &
&       (thermforces(:,:,:,inos)-qmass(inos)*dzeta(:,:,:,inos)*dzeta(:,:,:,inos+1))* &
&       dtion*dtion/qmass(inos)
   end if
 end do

 end select

 ABI_DEALLOCATE(thermforces)

end subroutine pimd_nosehoover_propagate
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_coord_transform
!! NAME
!!  pimd_coord_transform
!!
!! FUNCTION
!!  Apply a coordinate transformation on a given vector field
!!  (defined for each atom in in cell) - for PIMD
!!  Possible choices for the transformation:
!!    0: no transformation
!!    1: normal mode tranformation
!!    2: staging transformation
!!
!! INPUTS
!!  ioption=option given the direction of the transformation
!!     +1: from primitive coordinates to transformed coordinates
!!     -1: from transformed coordinates to primitive coordinates
!!  natom=number of atoms
!!  transform=coordinate transformation:
!!            0: no tranformation
!!            1: normal mode transformation
!!            2: staging transformation
!!  trotter=Trotter number
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  array(3,natom,trotter)=array to be transformed
!!
!! PARENTS
!!      pimd_langevin_nvt,pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_coord_transform(array,ioption,natom,transform,trotter)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ioption,natom,transform,trotter
!arrays
 real(dp),intent(inout) :: array(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,iimagem,iimagep,iimage2,ll
!arrays
 real(dp),allocatable :: array_temp(:,:,:)
 real(dp),allocatable :: nrm(:,:)

!************************************************************************

!=== No transformation ===================================================
 if (transform==0) then

   return

!=== Normal mode transformation ==========================================
 else if (transform==1) then

!  ---------------- From primitive to transformed coordinates ------------
   if (ioption==+1) then
     ABI_ALLOCATE(array_temp,(3,natom,trotter))

     array_temp(:,:,1)=zero
     do iimage=1,trotter
       array_temp(:,:,1)=array_temp(:,:,1)+array(:,:,iimage)
     end do

     do ll=2,trotter/2
       array_temp(:,:,2*ll-2)=zero
       array_temp(:,:,2*ll-1)=zero
       do iimage=1,trotter
         array_temp(:,:,2*ll-2)=array_temp(:,:,2*ll-2)+array(:,:,iimage)* &
&                              cos(two_pi*dble((iimage-1)*(ll-1))/dble(trotter))
         array_temp(:,:,2*ll-1)=array_temp(:,:,2*ll-1)-array(:,:,iimage)* &
&                              sin(two_pi*dble((iimage-1)*(ll-1))/dble(trotter))
       end do
     end do

     array_temp(:,:,trotter)=zero
     do iimage=1,trotter
       array_temp(:,:,trotter)=array_temp(:,:,trotter)+&
&                              array(:,:,iimage)*dble((-1)**(iimage-1))
     end do

     array=array_temp/dble(trotter)

     ABI_DEALLOCATE(array_temp)

!  ---------------- From transformed to primitive coordinates ------------
   else if (ioption==-1) then

    ABI_ALLOCATE(array_temp,(3,natom,trotter))  !real part
    ABI_ALLOCATE(nrm,(trotter,trotter))  !real part

   do iimage=1,trotter
     nrm(iimage,1)=one
     nrm(iimage,trotter)=dble((-1)**(iimage-1))
   end do

   do ll=2,trotter/2
     do iimage=1,trotter
       nrm(iimage,2*ll-2)= two*cos(two_pi*dble((ll-1)* &
&                        (iimage-1))/dble(trotter))
       nrm(iimage,2*ll-1)=-two*sin(two_pi*dble((ll-1)* &
&                        (iimage-1))/dble(trotter))
     end do
   end do

    do iimage=1,trotter
      array_temp(:,:,iimage)=zero
      do ll=1,trotter
        array_temp(:,:,iimage)=array_temp(:,:,iimage)+nrm(iimage,ll)*array(:,:,ll)
      end do
    end do
    array=array_temp

    ABI_DEALLOCATE(array_temp)
    ABI_DEALLOCATE(nrm)

   end if ! ioption

!=== Staging transformation ==============================================
 else if (transform==2) then

!  ---------------- From primitive to transformed coordinates ------------
   if (ioption==+1) then

     ABI_ALLOCATE(array_temp,(3,natom,trotter))
     array_temp=zero
     do iimage=1,trotter
       iimagep=iimage+1;if(iimage==trotter) iimagep=1
       iimagem=iimage-1;if(iimage==1) iimagem=trotter
       do iatom=1,natom
         do ii=1,3
           array_temp(ii,iatom,iimage)=(dble(iimagem)*array(ii,iatom,iimagep)+array(ii,iatom,1))/dble(iimage)
         end do
       end do
     end do
     if (trotter>1) then
       do iimage=2,trotter
         do iatom=1,natom
           do ii=1,3
             array(ii,iatom,iimage)=array(ii,iatom,iimage)-array_temp(ii,iatom,iimage)
           end do
         end do
       end do
     end if
     ABI_DEALLOCATE(array_temp)

!  ---------------- From transformed to primitive coordinates ------------
   else if (ioption==-1) then

     ABI_ALLOCATE(array_temp,(3,natom,trotter))
     array_temp=zero
     do iimage=1,trotter
       do iatom=1,natom
         do ii=1,3
           array_temp(ii,iatom,iimage)=array(ii,iatom,1)
         end do
       end do
     end do
     if (trotter>1) then
       do iimage=2,trotter
         do iatom=1,natom
           do ii=1,3
             do iimage2=iimage,trotter
               array_temp(ii,iatom,iimage)=array_temp(ii,iatom,iimage) &
&                     +array(ii,iatom,iimage2)*(dble(iimage-1))/(dble(iimage2-1))
             end do
           end do
         end do
       end do
     end if
     array(:,:,:)=array_temp(:,:,:)
     ABI_DEALLOCATE(array_temp)

   end if ! ioption

 end if ! transform

end subroutine pimd_coord_transform
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_force_transform
!! NAME
!!  pimd_force_transform
!!
!! FUNCTION
!!  Apply a coordinate transformation on forces (defined for each atom in in cell) - for PIMD
!!  Possible choices for the transformation:
!!    0: no transformation
!!    1: normal mode tranformation
!!    2: staging transformation
!!
!! INPUTS
!!  ioption=option given the direction of the transformation
!!     +1: from primitive coordinates to transformed coordinates
!!     -1: from transformed coordinates to primitive coordinates
!!  natom=number of atoms
!!  transform=coordinate transformation:
!!            0: no tranformation
!!            1: normal mode transformation
!!            2: staging transformation
!!  trotter=Trotter number
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  forces(3,natom,trotter)=array containing forces
!!
!! NOTES
!!  Back transformation (ioption=-1) not implemented !
!!
!! PARENTS
!!      pimd_langevin_nvt,pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_force_transform(forces,ioption,natom,transform,trotter)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ioption,natom,transform,trotter
!arrays
 real(dp),intent(inout) :: forces(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,ll
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: forces_temp(:,:,:),nrm(:,:)

!************************************************************************

 if (ioption==-1) then
   msg='Back transformation not implemented !'
   MSG_BUG(msg)
 end if

!=== No transformation ===================================================
 if (transform==0) then

   return

!=== Normal mode transformation ==========================================
 else if (transform==1) then

   !normal mode forces
   ABI_ALLOCATE(forces_temp,(3,natom,trotter))
   ABI_ALLOCATE(nrm,(trotter,trotter))

   do iimage=1,trotter
     nrm(iimage,1)=one
     nrm(iimage,trotter)=dble((-1)**(iimage-1))
   end do

   do ll=2,trotter/2
     do iimage=1,trotter
       nrm(iimage,2*ll-2)= two*cos(two_pi*dble((ll-1)* &
&                         (iimage-1))/dble(trotter))
       nrm(iimage,2*ll-1)=-two*sin(two_pi*dble((ll-1)* &
&                         (iimage-1))/dble(trotter))
     end do
   end do

   do ll=1,trotter
     forces_temp(:,:,ll)=zero
     do iimage=1,trotter
       forces_temp(:,:,ll)=forces_temp(:,:,ll)+nrm(iimage,ll)*forces(:,:,iimage)
     end do
   end do

   forces=forces_temp

   ABI_DEALLOCATE(forces_temp)
   ABI_DEALLOCATE(nrm)

!=== Staging transformation ==============================================
 else if (transform==2) then

   !staging forces
   ABI_ALLOCATE(forces_temp,(3,natom,trotter))
   forces_temp=zero
   do iimage=1,trotter
     do iatom=1,natom
       do ii=1,3
         forces_temp(ii,iatom,1)=forces_temp(ii,iatom,1)+forces(ii,iatom,iimage)
       end do
     end do
   end do
   if (trotter>1) then
     do iimage=2,trotter
       do iatom=1,natom
         do ii=1,3
           forces_temp(ii,iatom,iimage)=forces(ii,iatom,iimage) &
&               +forces_temp(ii,iatom,iimage-1)*(dble(iimage-2)/dble(iimage-1))
         end do
       end do
     end do
   end if
   forces=forces_temp
   ABI_DEALLOCATE(forces_temp)

 end if ! transform

end subroutine pimd_force_transform
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_apply_constraint
!! NAME
!!  pimd_apply_constraint
!!
!! FUNCTION
!!  Modify forces to take into account an holonomic constraint
!!  according to "pimd_constraint" parameter
!!  Available constraints:
!!    0: no constraint
!!    1: linear combination of coordinates
!!
!! INPUTS
!!  constraint=type of constraint to be applied
!!  mass(natom)=fictitious masses of atoms
!!  natom=number of atoms
!!  trotter=Trotter number
!!  wtatcon(3,natom)=weights for atomic constraints
!!  xcart(3,natom,trotter)=cartesian coordinates of atoms
!!
!! OUTPUT
!!  constraint_output(2)=several (real) data to be output
!!                       when a constraint has been applied
!!
!! SIDE EFFECTS
!!  forces(3,natom,trotter)=array containing forces
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt,pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_apply_constraint(constraint,constraint_output,forces,mass,natom,&
&                                trotter,wtatcon,xcart)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: constraint,natom,trotter
!arrays
 real(dp),intent(in) :: mass(natom),wtatcon(3,natom),xcart(3,natom,trotter)
 real(dp),intent(out) :: constraint_output(2)
 real(dp),intent(inout) :: forces(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage
 real(dp) :: af,lambda_cst,masstot,one_over_trotter,xcart_centroid,zz
 !character(len=500) :: msg
!arrays
 real(dp) :: force_centroid(3),lambda_com(3),mat(3,3),matinv(3,3),vec(3),weightsum(3)

!************************************************************************


 select case(constraint)

!=== No constraint =======================================================
 case(0)

   constraint_output(:)=zero
   return

!=== Linear combination of centroid coordinates ==========================
 case(1)

!  Some useful quantities
   masstot=sum(mass)*dble(trotter)
   weightsum(:)=sum(wtatcon,dim=2)*dble(trotter)
   zz=zero;af=zero;force_centroid=zero
   do iimage=1,trotter
     do iatom=1,natom
       do ii=1,3
         force_centroid(ii)=force_centroid(ii)+forces(ii,iatom,iimage)
         af=af+wtatcon(ii,iatom)*forces(ii,iatom,iimage)/mass(iatom)
         zz=zz+wtatcon(ii,iatom)**2/mass(iatom)
       end do
     end do
   end do
   vec(:)=force_centroid(:)-weightsum(:)*af/zz
   do ii=1,3
     mat(:,ii)=-weightsum(:)*weightsum(ii)/zz
     mat(ii,ii)=mat(ii,ii)+masstot
   end do
   call matr3inv(mat,matinv)

   !Calculation of a Lagrange multipliers:
   ! lambda_cst: to apply the constraint
   ! lambda_com: to maintain the position of the center of mass
   lambda_com(:)=matmul(matinv,vec)
   lambda_cst=(af-dot_product(weightsum,lambda_com))*dble(trotter)/zz

   !Modification of forces
   one_over_trotter=one/dble(trotter)
   do iimage=1,trotter
     do iatom=1,natom
       do ii=1,3
         forces(ii,iatom,iimage)=forces(ii,iatom,iimage) &
&                               -lambda_cst*wtatcon(ii,iatom)*one_over_trotter &
&                               -lambda_com(ii)*mass(iatom)
       end do
     end do
   end do

   !Computation of relevant outputs
   constraint_output(:)=zero
   !1-Reaction coordinate
   do iatom=1,natom
     do ii=1,3
       xcart_centroid=zero
       do iimage=1,trotter
         xcart_centroid=xcart_centroid+xcart(ii,iatom,iimage)
       end do
       constraint_output(1)=constraint_output(1)+xcart_centroid*wtatcon(ii,iatom)
     end do
   end do
   constraint_output(1)=constraint_output(1)/dble(trotter)
   !2-Force on reaction coordinate
   constraint_output(2)=-lambda_cst

 end select

end subroutine pimd_apply_constraint
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_mass_spring
!! NAME
!!  pimd_mass_spring
!!
!! FUNCTION
!!  Compute masses and spring constants for PIMD. Eventually apply a coordinate transformation.
!!  Possible choices for the transformation:
!!    0: no transformation
!!    1: normal mode tranformation
!!    2: staging transformation
!!
!! INPUTS
!!  inertmass(natom)=fictitious masses of atoms
!!  kt=kT constant
!!  natom=number of atoms
!!  quantummass(natom)=true masses of atoms
!!  transform=coordinate transformation:
!!            0: no tranformation
!!            1: normal mode transformation
!!            2: staging transformation
!!  trotter=Trotter number
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  spring(natom,mass_dim)=spring constants of atoms (mass_dim=1 or trotter)
!!
!! NOTES
!!  Back transformation (ioption=-1) not implemented !
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt,pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_mass_spring(inertmass,kt,mass,natom,quantummass,spring,transform,trotter)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,transform,trotter
 real(dp),intent(in) :: kt
!arrays
 real(dp),intent(in) :: inertmass(natom),quantummass(natom)
 real(dp),intent(out) :: mass(:,:),spring(:,:)
!Local variables-------------------------------
!scalars
 integer :: iimage,kk,natom_mass,natom_spring,nmass,nspring
 real(dp) :: gammasquare
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: mass_temp(:,:),lambda(:)

!************************************************************************

 natom_mass  =size(mass  ,1);nmass  =size(mass  ,2)
 natom_spring=size(spring,1);nspring=size(spring,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if
 if (natom/=natom_spring.or.(nspring/=1.and.nspring/=trotter)) then
   msg='Wrong dimensions for array spring !'
   MSG_BUG(msg)
 end if

!=== No transformation ===================================================
 if (transform==0) then
   !2nd dim of mass and spring = 1
   mass(1:natom,1)=inertmass(1:natom)
   spring(1:natom,1)=quantummass(1:natom)*dble(trotter)*kt*kt

!=== Normal mode transformation ==========================================
 else if (transform==1) then

   ABI_ALLOCATE(lambda,(trotter))
   lambda(1)=zero; lambda(trotter)=four*dble(trotter)
   do kk=2,trotter/2
     lambda(2*kk-2)=four*dble(trotter)* &
&                  (one-cos(two_pi*dble(kk-1)/dble(trotter)))
     lambda(2*kk-1)=four*dble(trotter)* &
&                  (one-cos(two_pi*dble(kk-1)/dble(trotter)))
   end do

   !normal mode masses
   do iimage=1,trotter
     mass(:,iimage)=quantummass(:)*lambda(iimage)
   end do

   do iimage=1,trotter
     spring(:,iimage)=mass(:,iimage)*dble(trotter)*kt*kt
   end do

   !fictitious masses
   mass(:,1)=inertmass(:)

   !from 2 to P not changed except if adiabatic PIMD
   !see Hone et al, JCP 124, 154103 (2006)
   gammasquare=one  !adiabaticity parameter
   do iimage=2,trotter
     mass(:,iimage)=mass(:,iimage)/gammasquare
   end do

   ABI_DEALLOCATE(lambda)

!=== Staging transformation ==============================================
 else if (transform==2) then

   !Fictitious masses
   mass(1:natom,1)=inertmass(1:natom)
   if (nmass>1) then
     do iimage=2,trotter
       mass(1:natom,iimage)=inertmass(1:natom)*dble(iimage)/dble(iimage-1)
     end do
   end if

   !Staging masses (mass_temp)
   ABI_ALLOCATE(mass_temp,(natom,trotter))
   mass_temp(1:natom,1)=zero
   if (nmass>1) then
     do iimage=2,trotter
       mass_temp(1:natom,iimage)=quantummass(1:natom)*dble(iimage)/dble(iimage-1)
     end do
   end if

   spring(1:natom,1)=mass_temp(1:natom,1)*dble(trotter)*kt*kt
   if (nspring>1) then
     do iimage=2,trotter
       spring(1:natom,iimage)=mass_temp(1:natom,iimage)*dble(trotter)*kt*kt
     end do
   end if
   ABI_DEALLOCATE(mass_temp)

 end if

end subroutine pimd_mass_spring
!!***

END MODULE m_pimd
!!***

!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_predtk
!! NAME
!!  m_predtk
!!
!! FUNCTION
!!  Low-level procedures used by 45_geomoptim routines
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, SE)
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

module m_predtk

 use defs_basis
 use m_profiling_abi
 use m_abimover

 implicit none

 private
!!***

 !public :: fdtion
!!***

contains
!!***

!!****f* ABINIT/fdtion
!! NAME
!! fdtion
!!
!! FUNCTION
!! Compute the apropiated "dtion" from the present values
!! of forces, velocity and viscosity
!!
!! INPUTS (in)
!! hist<type abihist>=Historical record of positions, forces
!!      |                    acell, stresses, and energies,
!! ab_mover<type abimover>=Subset of dtset only related with
!!          |                 movement of ions and acell, contains:
!!          | dtion:  Time step
!!          ! natom:  Number of atoms
!!          | vis:    viscosity
!!          | iatfix: Index of atoms and directions fixed
!!          | amass:  Mass of ions
!! itime: Index of time iteration
!! xcart(3,natom)= cartesian coordinates of atoms
!! fcart(3,natom)= forces in cartesian coordinates
!! vel(3,natom)= velocities
!!
!! OUTPUT (out)
!! fdtion = time step computed
!!
!! PARENTS
!!      pred_moldyn
!!
!! CHILDREN
!!
!! SOURCE

function fdtion(ab_mover,itime,xcart,fcart,vel)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fdtion'
!End of the abilint section

  implicit none

!Arguments ---------------------------------------------
!scalars
  type(abimover),intent(in) :: ab_mover
  integer,intent(in) :: itime
  real(dp) :: fdtion
!arrays
  real(dp) :: xcart(:,:),fcart(:,:),vel(:,:)

!Local variables ------------------------------
!scalars
  integer  :: jj,kk
  real(dp) :: max,min,val
  real(dp) :: ff,xc,vv,em

!************************************************************************

 max=0
 min=1e6

 do kk=1,ab_mover%natom
   em=ab_mover%amass(kk)
   do jj=1,3
     ff =fcart(jj,kk)
     xc =xcart(jj,kk)
     vv=vel(jj,kk)

     if (vv>1e-8) then
       val=abs(1.0_dp/vv)
       write(std_out,*) 'vel',kk,jj,val
       if (val>max) max=val
       if (val<min) min=val
     end if

     if (ff>1e-8) then
       val=sqrt(abs(2*em/ff))
       write(std_out,*) 'forces',kk,jj,val,em,ff
       if (val>max) max=val
       if (val<min) min=val
     end if

   end do

 end do

 write(std_out,*) "DTION max=",max
 write(std_out,*) "DTION min=",min

 if (itime==1)then
   fdtion=min/10
 else
   fdtion=min/10
 end if

 end function fdtion
!!***

end module m_predtk
!!***

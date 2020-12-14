!!****m* ABINIT/m_predtk
!! NAME
!!  m_predtk
!!
!! FUNCTION
!!  Low-level procedures used by 45_geomoptim routines
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, SE)
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
 use m_abicore
 use m_abimover

 implicit none

 private
!!***

 public :: fdtion
 public :: prtxvf        ! Print the values of xcart, vel, and fcart to unit iout.
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

!!****f* ABINIT/prtxvf
!!
!! NAME
!! prtxvf
!!
!! FUNCTION
!! Print the values of xcart, vel, and fcart to unit iout.
!! Also compute and print max and rms forces.

!! INPUTS
!! fcart(3,natom)=forces (hartree/bohr)
!! iatfix(3,natom)=1 for frozen or fixed atom along specified direction, else 0
!! iout=unit number for printing
!! natom=number of atoms in unit cell.
!! prtvel=1 to print velocities, else do not print them
!! vel(3,natom)=velocities
!! xcart(3,natom)=cartesian coordinates (bohr)
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      m_forces,m_gstateimg
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine prtxvf(fcart,fred,iatfix,iout,natom,prtvel,vel,xcart,xred)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,natom,prtvel
!arrays
 integer,intent(in) :: iatfix(3,natom)
 real(dp),intent(in) :: fcart(3,natom),fred(3,natom)
 real(dp),intent(in) :: xcart(3,natom),xred(3,natom)
 real(dp),intent(in) :: vel(3,natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,mu,unfixd
 real(dp) :: fmax,frms,val_max,val_rms
 character(len=500) :: msg

! *****************************************************************

 write(msg, '(a)' ) ' Cartesian coordinates (xcart) [bohr]'
 call wrtout(iout,msg,'COLL')
 do iatom=1,natom
   write(msg, '(1p,3e22.14)' )xcart(:,iatom)
   call wrtout(iout,msg,'COLL')
 end do

 write(msg, '(a)' ) ' Reduced coordinates (xred)'
 call wrtout(iout,msg,'COLL')
 do iatom=1,natom
   write(msg, '(1p,3e22.14)' )xred(:,iatom)
   call wrtout(iout,msg,'COLL')
 end do

!Compute max |f| and rms f, EXCLUDING the components determined by iatfix

 fmax=0.0_dp
 frms=0.0_dp
 unfixd=0
 do iatom=1,natom
   do mu=1,3
     if (iatfix(mu,iatom) /= 1) then
       unfixd=unfixd+1
       frms=frms+fcart(mu,iatom)**2
       fmax=max(fmax,abs(fcart(mu,iatom)))
     end if
   end do
 end do
 if ( unfixd /= 0 ) frms=sqrt(frms/dble(unfixd))

 write(msg, '(a,1p,2e12.5,a)' ) &
& ' Cartesian forces (fcart) [Ha/bohr]; max,rms=',fmax,frms,' (free atoms)'
 call wrtout(iout,msg,'COLL')
 do iatom=1,natom
   write(msg, '(1p,3e22.14)' )fcart(:,iatom)
   call wrtout(iout,msg,'COLL')
 end do

 write(msg, '(a)' ) ' Reduced forces (fred)'
 call wrtout(iout,msg,'COLL')
 do iatom=1,natom
   write(msg, '(1p,3e22.14)' )fred(:,iatom)
   call wrtout(iout,msg,'COLL')
 end do

 if (prtvel == 1) then

!  Compute max |v| and rms v,
!  EXCLUDING the components determined by iatfix
   val_max=0.0_dp
   val_rms=0.0_dp
   unfixd=0
   do iatom=1,natom
     do mu=1,3
       if (iatfix(mu,iatom) /= 1) then
         unfixd=unfixd+1
         val_rms=val_rms+vel(mu,iatom)**2
         val_max=max(val_max,abs(vel(mu,iatom)**2))
       end if
     end do
   end do
   if ( unfixd /= 0 ) val_rms=sqrt(val_rms/dble(unfixd))


   write(msg, '(a,1p,2e12.5,a)' ) ' Cartesian velocities (vel) [bohr*Ha/hbar]; max,rms=',&
&   sqrt(val_max),val_rms,' (free atoms)'
   call wrtout(iout,msg,'COLL')
   do iatom=1,natom
     write(msg, '(1p,3e22.14)' ) vel(:,iatom)
     call wrtout(iout,msg,'COLL')
   end do
 end if

end subroutine prtxvf
!!***

end module m_predtk
!!***

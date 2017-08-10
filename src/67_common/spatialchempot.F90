!{\src2tex{textfont=tt}}
!!****m* ABINIT/spatialchempot
!! NAME
!!  spatialchempot
!!
!! FUNCTION
!!  Treat spatially varying chemical potential.
!!  Compute energy and derivative with respect to dimensionless reduced atom coordinates of the
!!  spatially varying chemical potential. No contribution to stresses.
!!
!! COPYRIGHT
!!  Copyright (C) 2017 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! chempot(3,nzchempot,ntypat)=input array with information about the chemical potential (see input variable description)
!! natom=number of atoms in unit cell
!! ntypat=number of type of atoms
!! nzchempot=number of limiting planes for chemical potential
!! typat(natom)=integer label of each type of atom (1,2,...)
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!!
!! OUTPUT
!! e_chempot=chemical potential energy in hartrees
!! grchempottn(3,natom)=grads of e_chempot wrt xred(3,natom), hartrees.
!!
!! PARENTS
!!      setvtr
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine spatialchempot(e_chempot,chempot,grchempottn,natom,ntypat,nzchempot,typat,xred)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spatialchempot'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat,nzchempot
 real(dp),intent(out) :: e_chempot
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: chempot(3,nzchempot,ntypat),xred(3,natom)
 real(dp),intent(out) :: grchempottn(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,itypat,iz
 real(dp) :: a_2,a_3,cp0,cp1,dcp0,dcp1,ddz,deltaz,deltaziz
 real(dp) :: dqz,dz1,qz,zred,z0
!character(len=500) :: message

! *************************************************************************

!DEBUG
!write(std_out,'(a)')' spatialchempot : enter '
!write(std_out,'(a,2i6)')' nzchempot,ntypat',nzchempot,ntypat
!write(std_out,'(a,6es13.3)') ' chempot(1:3,1:2,1)=',chempot(1:3,1:2,1)
!write(std_out,'(a,6es13.3)') ' chempot(1:3,1:2,2)=',chempot(1:3,1:2,2)
!ENDDEBUG

 e_chempot=zero
 grchempottn(:,:)=zero

!Loop on the different atoms
 do iatom=1,natom

   itypat=typat(iatom)
   zred=xred(3,iatom)

!  Determine the delimiting plane just lower to zred
!  First compute the relative zred with respect to the first delimiting plane
!  Take into account a tolerance : 
   deltaz=zred-chempot(1,1,itypat)
   deltaz=modulo(deltaz+tol12,1.0d0)-tol12
!  deltaz is positive (or higher than -tol12), and lower than one-tol12.
   do iz=2,nzchempot+1
     if(iz/=nzchempot+1)then
       deltaziz=chempot(1,iz,itypat)-chempot(1,1,itypat)
     else
       deltaziz=one
     end if
     if(deltaziz>deltaz)exit
   end do

!  Defines coordinates and values inside the delimiting interval, 
!  with respect to the lower delimiting plane
   z0=chempot(1,iz-1,itypat)-chempot(1,1,itypat) ; cp0=chempot(2,iz-1,itypat) ; dcp0=chempot(3,iz-1,itypat)
   if(iz/=nzchempot+1)then
     dz1=chempot(1,iz,itypat)-chempot(1,iz-1,itypat) ; cp1=chempot(2,iz,itypat) ; dcp1=chempot(3,iz,itypat)
   else
     dz1=(chempot(1,1,itypat)+one)-chempot(1,nzchempot,itypat) ; cp1=chempot(2,1,itypat) ; dcp1=chempot(3,1,itypat)
   end if
   ddz=deltaz-z0

!DEBUG
!  write(std_out,'(a,2i5)')' Delimiting planes, iz-1 and iz=',iz-1,iz
!  write(std_out,'(a,2es13.3)')' z0,  dz1= :',z0,dz1
!  write(std_out,'(a,2es13.3)')' cp0, cp1= :',cp0,cp1
!  write(std_out,'(a,2es13.3)')' dcp0, dcp1= :',dcp0,dcp1
!  write(std_out,'(a,2es13.3)')' deltaz,ddz=',deltaz,ddz
!ENDDEBUG

!  Determine the coefficient of the third-order polynomial taking z0 as origin
!  P(dz=z-z0)= a_3*dz**3 + a_2*dz**2 + a_1*dz + a_0 ; obviously a_0=cp0 and a_1=dcp0
!  Define qz=a_3*dz + a_2 and dqz=3*a_3*dz + 2*a_2
   qz=((cp1-cp0)-dcp0*dz1)/dz1**2
   dqz=(dcp1-dcp0)/dz1
   a_3=(dqz-two*qz)/dz1
   a_2=three*qz-dqz

!  Compute value and gradient of the chemical potential, at ddz wrt to lower delimiting plane
   e_chempot=e_chempot+(a_3*ddz**3 + a_2*ddz**2 + dcp0*ddz + cp0)
   grchempottn(3,iatom)=three*a_3*ddz**2 + two*a_2*ddz + dcp0

!DEBUG
!  write(std_out,'(a,4es16.6)')' qz,dqz=',qz,dqz
!  write(std_out,'(a,4es16.6)')' cp0,dcp0,a_2,a_3=',cp0,dcp0,a_2,a_3
!  write(std_out,'(a,2es13.3)')' dcp0*ddz + cp0=',dcp0*ddz + cp0
!  write(std_out,'(a,2es13.3)')' a_2*ddz**2=',a_2*ddz**2
!  write(std_out,'(a,2es13.3)')' a_3*ddz**3=',a_3*ddz**3
!  write(std_out,'(a,2es13.3)')' contrib=',a_3*ddz**3 + a_2*ddz**2 + dcp0*ddz + cp0
!  write(std_out,'(a,2es13.3)')' e_chempot=',e_chempot
!  write(std_out,'(a,3es20.10)')' grchempottn=',grchempottn(:,iatom)
!ENDDEBUG

 end do

!DEBUG
!write(std_out,'(a)')' spatialchempot : exit '
!write(std_out,'(a,es16.6)') ' e_chempot=',e_chempot
!ENDDEBUG

end subroutine spatialchempot
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtxvf
!!
!! NAME
!! prtxvf
!!
!! FUNCTION
!! Print the values of xcart, vel, and fcart to unit iout.
!! Also compute and print max and rms forces.

!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, 
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! fcart(3,natom)=forces (hartree/bohr)
!! iatfix(3,natom)=1 for frozen or fixed atom along specified 
!!                 direction, else 0
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
!!      constrf,prtimg
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prtxvf(fcart,fred,iatfix,iout,natom,prtvel,vel,xcart,xred)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtxvf'
 use interfaces_14_hidewrite
!End of the abilint section

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
 character(len=500) :: message

! *****************************************************************

 write(message, '(a)' ) ' Cartesian coordinates (xcart) [bohr]'
 call wrtout(iout,message,'COLL')
 do iatom=1,natom
   write(message, '(1p,3e22.14)' )xcart(:,iatom)
   call wrtout(iout,message,'COLL')
 end do

 write(message, '(a)' ) ' Reduced coordinates (xred)'
 call wrtout(iout,message,'COLL')
 do iatom=1,natom
   write(message, '(1p,3e22.14)' )xred(:,iatom)
   call wrtout(iout,message,'COLL')
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

 write(message, '(a,1p,2e12.5,a)' ) &
& ' Cartesian forces (fcart) [Ha/bohr]; max,rms=',fmax,frms,' (free atoms)'
 call wrtout(iout,message,'COLL')
 do iatom=1,natom
   write(message, '(1p,3e22.14)' )fcart(:,iatom)
   call wrtout(iout,message,'COLL')
 end do

 write(message, '(a)' ) ' Reduced forces (fred)'
 call wrtout(iout,message,'COLL')
 do iatom=1,natom
   write(message, '(1p,3e22.14)' )fred(:,iatom)
   call wrtout(iout,message,'COLL')
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


   write(message, '(a,1p,2e12.5,a)' ) ' Cartesian velocities (vel) [bohr*Ha/hbar]; max,rms=',&
&   sqrt(val_max),val_rms,' (free atoms)'
   call wrtout(iout,message,'COLL')
   do iatom=1,natom
     write(message, '(1p,3e22.14)' ) vel(:,iatom)
     call wrtout(iout,message,'COLL')
   end do
 end if

end subroutine prtxvf
!!***

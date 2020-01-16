!!****m* ABINIT/m_pred_nose
!! NAME
!!  m_pred_nose
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR, JCC, SE)
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

module m_pred_nose

 use defs_basis
 use m_abicore
 use m_abimover
 use m_abihist

 use m_numeric_tools,  only : uniformrandom
 use m_geometry,    only : xcart2xred, xred2xcart, metric

 implicit none

 private
!!***

 public :: pred_nose
!!***

contains
!!***

!!****f* ABINIT/pred_nose
!! NAME
!! pred_nose
!!
!! FUNCTION
!! Ionmov predictors (8) Verlet algorithm with a nose-hoover thermostat
!!
!! IONMOV 8:
!! Given a starting point xred that is a vector of length 3*natom
!! (reduced nuclei coordinates), a velocity vector (in cartesian
!! coordinates), and unit cell parameters (acell and rprimd -
!! without velocities in the present implementation),
!! the Verlet dynamics is performed, using the gradient of the
!! energy (atomic forces and stresses) as calculated by the routine scfcv.
!!
!! Some atoms can be kept fixed, while the propagation of unit cell
!! parameters is only performed if optcell/=0.
!! No more than "ntime" steps are performed.
!! The time step is governed by dtion (contained in dtset)
!! Returned quantities are xred, and eventually acell and rprimd
!! (new ones!).
!!
!! See ionmov=6, but with a nose-hoover thermostat
!! Velocity verlet algorithm : Swope et al JCP 76 (1982) 637
!!
!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information
!!                                needed by the preditor
!! itime  : Index of the present iteration
!! ntime  : Maximal number of iterations
!! zDEBUG : if true print some debugging information
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces acell, rprimd, stresses
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      hist2var,metric,var2hist,wrtout,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine pred_nose(ab_mover,hist,itime,ntime,zDEBUG,iexit)

 implicit none

!Arguments ------------------------------------
!scalars
 type(abimover),intent(in)       :: ab_mover
 type(abihist),intent(inout) :: hist
 integer,intent(in) :: itime
 integer,intent(in) :: ntime
 integer,intent(in) :: iexit
 logical,intent(in) :: zDEBUG

!Local variables-------------------------------
!scalars
 integer  :: ii,jj,kk
 integer  :: idum=-5
 real(dp),parameter :: v2tol=tol8,nosetol=tol10
 real(dp) :: delxi,xio,ktemp,rescale_vel
 real(dp) :: dnose,v2nose,xin_nose
 real(dp),save :: xi_nose,fsnose,snose
 real(dp) :: gnose
 real(dp) :: ucvol,ucvol_next
 real(dp) :: etotal
 logical  :: ready

!arrays
 real(dp) :: acell(3),acell_next(3)
 real(dp) :: rprimd(3,3),rprimd_next(3,3)
 real(dp) :: gprimd(3,3)
 real(dp) :: gmet(3,3)
 real(dp) :: rmet(3,3)
 real(dp) :: fcart(3,ab_mover%natom)
!real(dp) :: fred_corrected(3,ab_mover%natom)
 real(dp) :: xred(3,ab_mover%natom),xred_next(3,ab_mover%natom)
 real(dp) :: xcart(3,ab_mover%natom),xcart_next(3,ab_mover%natom)
 real(dp) :: vel(3,ab_mover%natom),vel_temp(3,ab_mover%natom)
 real(dp) :: finose(3,ab_mover%natom),binose(3,ab_mover%natom)
 real(dp) :: vonose(3,ab_mover%natom),hnose(3,ab_mover%natom)
 real(dp),allocatable,save :: fcart_m(:,:),fcart_mold(:,:)
 real(dp) :: strten(6)
 character(len=500) :: message

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 if(iexit/=0)then
   if (allocated(fcart_m))     then
     ABI_DEALLOCATE(fcart_m)
   end if
   if (allocated(fcart_mold))  then
     ABI_DEALLOCATE(fcart_mold)
   end if
   return
 end if

!write(std_out,*) 'nose 01'
!##########################################################
!### 01. Allocate the arrays fcart_m and fcart_mold

 if(itime==1)then
   if (allocated(fcart_m))     then
     ABI_DEALLOCATE(fcart_m)
   end if
   if (allocated(fcart_mold))  then
     ABI_DEALLOCATE(fcart_mold)
   end if
 end if

 if(.not.allocated(fcart_m))     then
   ABI_ALLOCATE(fcart_m,(3,ab_mover%natom))
 end if
 if(.not.allocated(fcart_mold))  then
   ABI_ALLOCATE(fcart_mold,(3,ab_mover%natom))
 end if

!write(std_out,*) 'nose 02'
!##########################################################
!### 02. Obtain the present values from the history

 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 fcart(:,:)=hist%fcart(:,:,hist%ihist)
 strten(:)=hist%strten(:,hist%ihist)
 vel(:,:)=hist%vel(:,:,hist%ihist)
 etotal=hist%etot(hist%ihist)

 write(std_out,*) 'RPRIMD'
 do ii=1,3
   write(std_out,*) rprimd(:,ii)
 end do
 write(std_out,*) 'RMET'
 do ii=1,3
   write(std_out,*) rmet(ii,:)
 end do

!Get rid of mean force on whole unit cell, but only if no
!generalized constraints are in effect
!  call fcart2fred(fcart,fred_corrected,rprimd,ab_mover%natom)
!  if(ab_mover%nconeq==0)then
!    amass_tot=sum(ab_mover%amass(:))
!    do ii=1,3
!      if (ii/=3.or.ab_mover%jellslab==0) then
!        favg=sum(fred_corrected(ii,:))/dble(ab_mover%natom)
!        fred_corrected(ii,:)=fred_corrected(ii,:)-favg*ab_mover%amass(:)/amass_tot
!      end if
!    end do
!  end if

!write(std_out,*) 'nose 03'
!##########################################################
!### 03. Fill the vectors vin and vout

!write(std_out,*) 'nose 04'
!##########################################################
!### 04. Initialize or update the hessian matrix

!write(std_out,*) 'nose 05'
!##########################################################
!### 05. Compute the next values

!The temperature is linear between initial and final values
!It is here converted from Kelvin to Hartree (kb_HaK)
 ktemp=(ab_mover%mdtemp(1)+((ab_mover%mdtemp(2)-ab_mover%mdtemp(1))/dble(ntime-1))*(itime-1))*kb_HaK

!%%% NOSE DYNAMICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 acell_next(:)=acell(:)
 ucvol_next=ucvol
 rprimd_next(:,:)=rprimd(:,:)

 if(itime==1)then
   snose=0.0_dp
   xi_nose=0.0_dp
!  Compute twice the kinetic energy of the system, called v2nose
   v2nose=0.0_dp
   do kk=1,ab_mover%natom
     do jj=1,3
       v2nose=v2nose+vel(jj,kk)*vel(jj,kk)*ab_mover%amass(kk)
     end do
   end do
   if (zDEBUG)then
     write(std_out,*) 'itime ntime KTEMP=',itime-1,ntime-1,ktemp
     write(std_out,*) 'V2NOSE=',v2nose
     write (std_out,*) 'VEL'
     do kk=1,ab_mover%natom
       write (std_out,*) vel(:,kk)
     end do
   end if

!  If there is no kinetic energy, use a random initial velocity
   if (v2nose<=v2tol) then
     v2nose=0.0_dp
     do kk=1,ab_mover%natom
       do jj=1,3
!        Uniform random returns a uniform random deviate between 0.0
!        and 1.0
!        if it were always 0 or 1, then the following expression
!        would give the requested temperature
         vel(jj,kk)=(1.0_dp-2.0_dp*uniformrandom(idum))*&
&         sqrt( (ab_mover%mdtemp(1)) * kb_HaK / ab_mover%amass(kk) )
!        Recompute v2nose
         v2nose=v2nose+vel(jj,kk)*vel(jj,kk)*ab_mover%amass(kk)
         if (zDEBUG)then
           write(std_out,*) 'jj kk vel(jj,kk)=',jj,kk,vel(jj,kk)
           write(std_out,*) 'jj kk V2NOSE=',jj,kk,v2nose
         end if
       end do
     end do
   end if
   write(std_out,*) 'V2NOSE=',v2nose

!  Now, rescale the velocities to give the proper temperature
   rescale_vel=sqrt(3.0_dp*ab_mover%natom*(ab_mover%mdtemp(1))*kb_HaK/v2nose)
   write(std_out,*) 'RESCALE_VEL=',rescale_vel
   vel(:,:)=vel(:,:)*rescale_vel
!  Recompute v2nose with the rescaled velocities
   v2nose=0.0_dp
   do kk=1,ab_mover%natom
     do jj=1,3
       v2nose=v2nose+vel(jj,kk)*vel(jj,kk)*ab_mover%amass(kk)
     end do
   end do
   write(message, '(a)' )&
&   ' Rescaling or initializing velocities to initial temperature'
   call wrtout(std_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(2(a,es22.14))' )&
&   ' ---  Scaling factor : ',rescale_vel,&
&   ' Asked T (K) ',ab_mover%mdtemp(1)
   call wrtout(std_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es22.14)' )&
&   ' ---  Effective temperature',v2nose/(3.0_dp*ab_mover%natom*kb_HaK)
   call wrtout(std_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 do kk=1,ab_mover%natom
   do jj=1,3
     fcart_m(jj,kk)=fcart(jj,kk)/ab_mover%amass(kk)
   end do
 end do

!First step of velocity verlet algorithm
 gnose=3*ab_mover%natom

!Convert input xred (reduced coordinates) to xcart (cartesian)
 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

!Calculate nose-hoover force on atoms
!If first iteration, no old force are available, so use present
!forces
 if (itime==1) fcart_mold(:,:)=fcart_m(:,:)

 if (zDEBUG)then
   write (std_out,*) 'FCART_MOLD'
   do kk=1,ab_mover%natom
     write (std_out,*) fcart_mold(:,kk)
   end do
   write (std_out,*) 'FCART_M'
   do kk=1,ab_mover%natom
     write (std_out,*) fcart_m(:,kk)
   end do
 end if

 finose(:,:)=fcart_mold(:,:)-xi_nose*vel(:,:)
 xcart(:,:)=xcart(:,:)+ab_mover%dtion*(vel(:,:)+ab_mover%dtion*finose(:,:)/2.0_dp)

!Convert back to xred (reduced coordinates)
 call xcart2xred(ab_mover%natom,rprimd,xcart,xred)

 if (zDEBUG)then
   write (std_out,*) 'VEL'
   do kk=1,ab_mover%natom
     write (std_out,*) vel(:,kk)
   end do
 end if

!Calculate v2nose
 v2nose=0.0_dp
 do kk=1,ab_mover%natom
   do jj=1,3
     v2nose=v2nose+vel(jj,kk)*vel(jj,kk)*ab_mover%amass(kk)
   end do
 end do
 vel(:,:)=vel(:,:)+ab_mover%dtion*finose(:,:)/2.0_dp

 if (zDEBUG)then
   write(std_out,*) 'NOSE BEFORE'
   write(std_out,*) 'FSNOSE=',fsnose
   write(std_out,*) 'SNOSE=',snose
   write(std_out,*) 'XI_NOSE=',xi_nose
   write (std_out,*) 'VEL'
   do kk=1,ab_mover%natom
     write (std_out,*) vel(:,kk)
   end do
   write (std_out,*) 'NOSEINERT',ab_mover%noseinert
 end if

!Update thermostat
 fsnose=(v2nose-gnose*ktemp)/ab_mover%noseinert
 snose=snose+ab_mover%dtion*(xi_nose+ab_mover%dtion*fsnose/2.0_dp)
 xi_nose=xi_nose+ab_mover%dtion*fsnose/2.0_dp
 if (zDEBUG)then
   write(std_out,*) 'NOSE AFTER'
   write(std_out,*) 'FSNOSE=',fsnose
   write(std_out,*) 'SNOSE=',snose
   write(std_out,*) 'XI_NOSE=',xi_nose
   write (std_out,*) 'VEL'
   do kk=1,ab_mover%natom
     write (std_out,*) vel(:,kk)
   end do
 end if

!Second step of the velocity Verlet algorithm, uses the 'new forces'
!Calculate v2nose
 v2nose=0.0_dp
 do kk=1,ab_mover%natom
   do jj=1,3
     v2nose=v2nose+vel(jj,kk)*vel(jj,kk)*ab_mover%amass(kk)
   end do
 end do
 vel_temp(:,:)=vel(:,:)

 if (zDEBUG)then
   write(std_out,*) 'V2NOSE=',v2nose
   write (std_out,*) 'VEL'
   do kk=1,ab_mover%natom
     write (std_out,*) vel(:,kk)
   end do
   write (std_out,*) 'Starting Newton Raphson'
 end if

 xin_nose=xi_nose

!Start Newton-Raphson loop
 ready=.false.
 do while (.not.ready)
   xio=xin_nose
   delxi=0.0D0
   vonose(:,:)=vel_temp(:,:)
   hnose(:,:)=-ab_mover%dtion/2.0_dp*(fcart_m(:,:)-xio*vonose(:,:))-&
&   (vel(:,:)-vonose(:,:))
   do kk=1,ab_mover%natom
     do jj=1,3
       binose(jj,kk)=vonose(jj,kk)*ab_mover%dtion/ab_mover%noseinert*ab_mover%amass(kk) ! a verifier
       delxi=delxi+hnose(jj,kk)*binose(jj,kk)
     end do
   end do
   dnose=-(xio*ab_mover%dtion/2.0D0+1.0D0)
   delxi=delxi-dnose*((-v2nose+gnose*ktemp)*ab_mover%dtion/2.0_dp/ &
&   ab_mover%noseinert-(xi_nose-xio))
   delxi=delxi/(-ab_mover%dtion*ab_mover%dtion/2.0_dp*v2nose/ab_mover%noseinert+dnose)

!  hzeronose=-(xio-xi_nose-(v2nose-gnose*ktemp)
!  *dtion/(2.0_dp*ab_mover%noseinert) )
!  cibinose=-v2nose*dtion*dtion/(2.0_dp*ab_mover%noseinert)
!  delxi=(delxi+hzeronose*dnose)/(dnose+cibinose)

!  DEBUG
!  write(message, '(a,es22.14)' )' after delxi',delxi
!  call wrtout(std_out,message,'COLL')
!  call wrtout(std_out,message,'COLL')
!  ENDDEBUG
   v2nose=0.0_dp

   vel_temp(:,:)=vel_temp(:,:)+&
&   (hnose+ab_mover%dtion/2.0_dp*vonose(:,:)*delxi)/dnose
   do kk=1,ab_mover%natom
     do jj=1,3
       v2nose=v2nose+vel_temp(jj,kk)*&
&       vel_temp(jj,kk)*ab_mover%amass(kk)
     end do
   end do
!  New guess for xi
   xin_nose=xio+delxi

!  zDEBUG
!  write(message, '(a,es22.14)' )' v2nose=',v2nose
!  call wrtout(std_out,message,'COLL')
!  call wrtout(std_out,message,'COLL')
!  ENDDEBUG

   ready=.true.
!  Test for convergence
   kk=0
   jj=1
   do while((kk<=ab_mover%natom).and.(jj<=3).and.ready)
     kk=kk+1
     if (kk>ab_mover%natom) then
       kk=1
       jj=jj+1
     end if
     if ((kk<=ab_mover%natom) .and.(jj<=3)) then
       if (abs(vel_temp(jj,kk))<1.0d-50)&
&       vel_temp(jj,kk)=1.0d-50
       if (abs((vel_temp(jj,kk)-vonose(jj,kk))&
&       /vel_temp(jj,kk))>nosetol) ready=.false.
     else
       if (xin_nose<1.0d-50) xin_nose=1.0d-50
       if (abs((xin_nose-xio)/xin_nose)>nosetol) ready=.false.
     end if
   end do   ! end of while

!  Enddo ready
 end do

!Update velocities to converged value
 vel(:,:)=vel_temp(:,:)
 write(message, '(a,es14.7)' )' converged velocities for T=',ktemp
 call wrtout(std_out,message,'COLL')

 if (zDEBUG)then
   write (std_out,*) 'Final Values for NOSE'
   write (std_out,*) 'VEL'
   do kk=1,ab_mover%natom
     write (std_out,*) vel(:,kk)
   end do
   write (std_out,*) 'XCART'
   do kk=1,ab_mover%natom
     write (std_out,*) xcart(:,kk)
   end do
 end if

!Update thermostat
 xi_nose=xin_nose
 xcart_next(:,:)=xcart(:,:)
!Convert back to xred_next (reduced coordinates)
 call xcart2xred(ab_mover%natom,rprimd,xcart_next,xred_next)
!Store 'new force' as 'old force'
 fcart_mold(:,:)=fcart_m(:,:)

!write(std_out,*) 'nose 06'
!##########################################################
!### 06. Update the history with the prediction

!Increase indexes
 hist%ihist=abihist_findIndex(hist,+1)

 call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 hist%vel(:,:,hist%ihist)=vel(:,:)
 hist%time(hist%ihist)=real(itime,kind=dp)*ab_mover%dtion

end subroutine pred_nose
!!***

end module m_pred_nose
!!***

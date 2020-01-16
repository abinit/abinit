!!****m* ABINIT/m_pred_isokinetic
!! NAME
!!  m_pred_isokinetic
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

module m_pred_isokinetic

 use m_abicore
 use defs_basis
 use m_abimover
 use m_abihist

 use m_numeric_tools,  only : uniformrandom
 use m_geometry,  only : xcart2xred, xred2xcart


 implicit none

 private
!!***

 public :: pred_isokinetic
!!***

contains
!!***


!!****f* ABINIT/pred_isokinetic
!! NAME
!! pred_isokinetic
!!
!! FUNCTION
!! Ionmov predictors (12) Isokinetic ensemble molecular dynamics
!!
!! IONMOV 12:
!! Isokinetic ensemble molecular dynamics.
!! The equation of motion of the ions in contact with a thermostat
!! are solved with the algorithm proposed by Zhang [J. Chem. Phys. 106, 6102 (1997)] [[cite:Zhang1997]],
!! as worked out by Minary et al, J. Chem. Phys. 188, 2510 (2003) [[cite:Minary2003]].
!! The conservation of the kinetic energy is obtained within machine precision, at each step.
!! Related parameters : the time step (dtion), the initial temperature (mdtemp(1)) if the velocities are not defined to start with.
!!
!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information needed by the preditor
!! itime  : Index of the present iteration
!! ntime  : Maximal number of iterations
!! zDEBUG : if true print some debugging information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces acell, rprimd, stresses
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      hist2var,var2hist,wrtout,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine pred_isokinetic(ab_mover,hist,itime,ntime,zDEBUG,iexit)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itime
 integer,intent(in) :: ntime
 integer,intent(in) :: iexit
 logical,intent(in) :: zDEBUG
 type(abimover),intent(in)       :: ab_mover
 type(abihist),intent(inout) :: hist

!Local variables-------------------------------
!scalars
 integer  :: kk,iatom,idim,idum=5,nxyzatfree,ndegfreedom,nfirst,ifirst
 real(dp) :: a,as,b,sqb,s,s1,s2,scdot,sigma2,vtest,v2gauss
 real(dp),parameter :: v2tol=tol8
 real(dp) :: etotal,rescale_vel
 character(len=5000) :: message
!arrays
 real(dp),allocatable,save :: fcart_m(:,:),vel_nexthalf(:,:)

 real(dp) :: acell(3),rprimd(3,3)
 real(dp) :: fcart(3,ab_mover%natom)
!real(dp) :: fred_corrected(3,ab_mover%natom)
 real(dp) :: xcart(3,ab_mover%natom),xcart_next(3,ab_mover%natom)
 real(dp) :: xred(3,ab_mover%natom),xred_next(3,ab_mover%natom)
 real(dp) :: vel(3,ab_mover%natom)
 real(dp) :: strten(6)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

!DEBUG
!write(std_out,*)' pred_isokinetic : enter '
!stop
!ENDDEBUG

 if(iexit/=0)then
   if (allocated(fcart_m))       then
     ABI_DEALLOCATE(fcart_m)
   end if
   if (allocated(vel_nexthalf))  then
     ABI_DEALLOCATE(vel_nexthalf)
   end if
   return
 end if

!write(std_out,*) 'isokinetic 01'
!##########################################################
!### 01. Debugging and Verbose

 if(zDEBUG)then
   write(std_out,'(a,3a,40a,37a)') ch10,('-',kk=1,3),&
&   'Debugging and Verbose for pred_isokinetic',('-',kk=1,37)
   write(std_out,*) 'ionmov: ',12
   write(std_out,*) 'itime:  ',itime
 end if

!write(std_out,*) 'isokinetic 02'
!##########################################################
!### 02. Allocate the vectors vin, vout and hessian matrix
!###     These arrays could be allocated from a previous
!###     dataset that exit before itime==ntime

 if(itime==1)then
   if (allocated(fcart_m))       then
     ABI_DEALLOCATE(fcart_m)
   end if
   if (allocated(vel_nexthalf))  then
     ABI_DEALLOCATE(vel_nexthalf)
   end if
 end if

 if (.not.allocated(fcart_m))       then
   ABI_ALLOCATE(fcart_m,(3,ab_mover%natom))
 end if
 if (.not.allocated(vel_nexthalf))  then
   ABI_ALLOCATE(vel_nexthalf,(3,ab_mover%natom))
 end if

!write(std_out,*) 'isokinetic 03'
!##########################################################
!### 03. Obtain the present values from the history

 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)

 fcart(:,:)=hist%fcart(:,:,hist%ihist)
 strten(:) =hist%strten(:,hist%ihist)
 vel(:,:)  =hist%vel(:,:,hist%ihist)
 etotal    =hist%etot(hist%ihist)

 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

 if(zDEBUG)then
   write (std_out,*) 'fcart:'
   do kk=1,ab_mover%natom
     write (std_out,*) fcart(:,kk)
   end do
   write (std_out,*) 'vel:'
   do kk=1,ab_mover%natom
     write (std_out,*) vel(:,kk)
   end do
   write (std_out,*) 'strten:'
   write (std_out,*) strten(1:3),ch10,strten(4:6)
   write (std_out,*) 'etotal:'
   write (std_out,*) etotal
 end if

!Get rid of mean force on whole unit cell, but only if no
!generalized constraints are in effect
!  call fcart2fred(hist%fcart(:,:,hist%ihist),fred_corrected,rprimd,ab_mover%natom)
!  if(ab_mover%nconeq==0)then
!    amass_tot=sum(ab_mover%amass(:))
!    do ii=1,3
!      if (ii/=3.or.ab_mover%jellslab==0) then
!        favg=sum(fred_corrected(ii,:))/dble(ab_mover%natom)
!        fred_corrected(ii,:)=fred_corrected(ii,:)-favg*ab_mover%amass(:)/amass_tot
!      end if
!    end do
!  end if

!Count the number of degrees of freedom, taking into account iatfix.
!Also fix the velocity to zero for the fixed atoms
 nxyzatfree=0
 do iatom=1,ab_mover%natom
   do idim=1,3
     if(ab_mover%iatfix(idim,iatom)==0)then
       nxyzatfree=nxyzatfree+1
     else
       vel(idim,iatom)=zero
     endif
   enddo
 enddo

!Now, the number of degrees of freedom is reduced by four because of the kinetic energy conservation
!and because of the conservation of the total momentum for each dimension, in case no atom position is fixed for that dimension 
!(in the latter case, one degree of freedom has already been taken away)
!This was not done until v8.9 of ABINIT ...
 ndegfreedom=nxyzatfree
 ndegfreedom=nxyzatfree-1 ! Kinetic energy conservation
 do idim=1,3
   if(sum(ab_mover%iatfix(idim,:))==0)then
     ndegfreedom=ndegfreedom-1 ! Macroscopic momentum
   endif
 enddo

!write(std_out,*) 'isokinetic 04'
!##########################################################
!### 04. Second half-velocity (Only after the first itime)

 if(itime>1) then

   do iatom=1,ab_mover%natom
     do idim=1,3
       if(ab_mover%iatfix(idim,iatom)==0)then
         fcart_m(idim,iatom)=fcart(idim,iatom)/ab_mover%amass(iatom)
       else
         fcart_m(idim,iatom)=zero
       endif
     end do
   end do

!  Computation of vel(:,:) at the next positions
!  Computation of v2gauss, actually twice the kinetic energy. 
!  Called 2K, cf Eq. (A13) of [[cite:Minary2003]].
   v2gauss=0.0_dp
   do iatom=1,ab_mover%natom
     do idim=1,3
       v2gauss=v2gauss+&
&       vel_nexthalf(idim,iatom)*vel_nexthalf(idim,iatom)*&
&       ab_mover%amass(iatom)
     end do
   end do

!  Computation of a and b (4.13 of [[cite:Minary2003]])
   a=0.0_dp
   b=0.0_dp
   do iatom=1,ab_mover%natom
     do idim=1,3
       a=a+fcart_m(idim,iatom)*vel_nexthalf(idim,iatom)*ab_mover%amass(iatom)
       b=b+fcart_m(idim,iatom)*fcart_m(idim,iatom)*ab_mover%amass(iatom)
     end do
   end do
   a=a/v2gauss
   b=b/v2gauss


!  Computation of s and scdot
   sqb=sqrt(b)
   as=sqb*ab_mover%dtion/2.
! jmb
   if ( as > 300.0 ) as=300.0
   s1=cosh(as)
   s2=sinh(as)
   s=a*(s1-1.)/b+s2/sqb
   scdot=a*s2/sqb+s1

   do iatom=1,ab_mover%natom
     do idim=1,3
       if(ab_mover%iatfix(idim,iatom)==0)then
         vel(idim,iatom)=(vel_nexthalf(idim,iatom)+fcart_m(idim,iatom)*s)/scdot
       else
         vel(idim,iatom)=zero
       endif
     enddo
   enddo  

   if (zDEBUG)then
     write(std_out,*) 'Computation of the second half-velocity'
     write(std_out,*) 'Cartesian forces per atomic mass (fcart_m):'
     do kk=1,ab_mover%natom
       write (std_out,*) fcart_m(:,kk)
     end do
     write(std_out,*) 'vel:'
     do kk=1,ab_mover%natom
       write (std_out,*) vel(:,kk)
     end do
     write(std_out,*) 'v2gauss:',v2gauss
     write(std_out,*) 'a:',a
     write(std_out,*) 'b:',b
     write(std_out,*) 's:',s
     write(std_out,*) 'scdot:',scdot
   end if

 end if ! (if itime>1)

!write(std_out,*) 'isokinetic 05'
!##########################################################
!### 05. First half-time (First cycle the loop is double)

 if (itime==1) then
   nfirst=2
 else
   nfirst=1
 end if

 do ifirst=1,nfirst

!  Application of Gauss' principle of least constraint according to Fei Zhang's algorithm (J. Chem. Phys. 106, 1997, [[cite:Zhang1997]] p.6102)

!  v2gauss is twice the kinetic energy
   v2gauss=0.0_dp
   do iatom=1,ab_mover%natom
     do idim=1,3
       v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*ab_mover%amass(iatom)
     end do
   end do

!  If there is no kinetic energy to start with ...
   if (v2gauss<=v2tol.and.itime==1) then
!    Maxwell-Boltzman distribution
     v2gauss=zero
     vtest=zero
     do iatom=1,ab_mover%natom
       do idim=1,3
         if(ab_mover%iatfix(idim,iatom)==0)then
           vel(idim,iatom)=sqrt(kb_HaK*ab_mover%mdtemp(1)/ab_mover%amass(iatom))*cos(two_pi*uniformrandom(idum))
           vel(idim,iatom)=vel(idim,iatom)*sqrt(-2._dp*log(uniformrandom(idum)))
         else
           vel(idim,iatom)=zero
         endif
       end do
     end do

!    Get rid of center-of-mass velocity
     s1=sum(ab_mover%amass(:))
     do idim=1,3
       if(sum(ab_mover%iatfix(idim,:))==0)then
         s2=sum(ab_mover%amass(:)*vel(idim,:))
         vel(idim,:)=vel(idim,:)-s2/s1
       endif
     end do

!    Recompute v2gauss
     do iatom=1,ab_mover%natom
       do idim=1,3
         v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*ab_mover%amass(iatom)
         vtest=vtest+vel(idim,iatom)/ndegfreedom
       end do
     end do

!    Now rescale the velocities to give the exact temperature
     rescale_vel=sqrt(ndegfreedom*kb_HaK*ab_mover%mdtemp(1)/v2gauss)
     vel(:,:)=vel(:,:)*rescale_vel

!    Recompute v2gauss with the rescaled velocities
     v2gauss=zero
     do iatom=1,ab_mover%natom
       do idim=1,3
         v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*ab_mover%amass(iatom)
       end do
     end do

!    Compute the variance and print
     sigma2=(v2gauss/ndegfreedom-ab_mover%amass(1)*vtest**2)/kb_HaK

   end if

   do iatom=1,ab_mover%natom
     do idim=1,3
       if(ab_mover%iatfix(idim,iatom)==0)then
         fcart_m(idim,iatom)=fcart(idim,iatom)/ab_mover%amass(iatom)
       else
         fcart_m(idim,iatom)=zero
       endif
     end do
   end do

   if (zDEBUG)then
     write(std_out,*) 'Calculation first half-velocity '
     write (std_out,*) 'vel:'
     do kk=1,ab_mover%natom
       write (std_out,*) vel(:,kk)
     end do
     write (std_out,*) 'xcart:'
     do kk=1,ab_mover%natom
       write (std_out,*) xcart(:,kk)
     end do
     write (std_out,*) 'xred:'
     do kk=1,ab_mover%natom
       write (std_out,*) xred(:,kk)
     end do
     write (std_out,*) 'fcart_m'
     do kk=1,ab_mover%natom
       write (std_out,*) fcart_m(:,kk)
     end do
     write(std_out,*) 's2',s2
     write(std_out,*) 'v2gauss',v2gauss
     write(std_out,*) 'sigma2',sigma2

     write(message, '(a)' )&
&     ' --- Rescaling or initializing velocities to initial temperature'
     call wrtout(std_out,message,'COLL')
     write(message, '(a,d12.5,a,D12.5)' )&
&     ' --- Scaling factor :',rescale_vel,' Asked T (K) ',ab_mover%mdtemp(1)
     call wrtout(std_out,message,'COLL')
     write(message, '(a,d12.5,a,D12.5)' )&
&     ' --- Effective temperature',v2gauss/(ndegfreedom*kb_HaK),' From variance', sigma2
     call wrtout(std_out,message,'COLL')
   end if

!  Convert input xred (reduced coordinates) to xcart (cartesian)
   call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

   if(itime==1.and.ifirst==1) then
     call wrtout(std_out,'if itime==1','COLL')
     vel_nexthalf(:,:)=vel(:,:)
     xcart_next(:,:)=xcart(:,:)
     call xcart2xred(ab_mover%natom,rprimd,xcart_next,xred_next)
     xred=xred_next
     call xred2xcart(ab_mover%natom,rprimd,xcart,xred)
   end if

 end do

!Computation of vel_nexthalf (4.16 of [[cite:Minary2003]])
!Computation of a and b (4.13 of [[cite:Minary2003]])
 a=0.0_dp
 b=0.0_dp
 do iatom=1,ab_mover%natom
   do idim=1,3
     a=a+fcart_m(idim,iatom)*vel(idim,iatom)*ab_mover%amass(iatom)
     b=b+fcart_m(idim,iatom)*fcart_m(idim,iatom)*ab_mover%amass(iatom)
   end do
 end do
 a=a/v2gauss+tol20
 b=b/v2gauss+tol20
!Computation of s and scdot
 sqb=sqrt(b)+tol20
 as=sqb*ab_mover%dtion/2.
! jmb
 if ( as > 300.0 ) as=300.0
 s1=cosh(as)
 s2=sinh(as)
 s=a*(s1-1.)/b+s2/sqb
 scdot=a*s2/sqb+s1
 do iatom=1,ab_mover%natom
   do idim=1,3
     if(ab_mover%iatfix(idim,iatom)==0)then
       vel_nexthalf(idim,iatom)=(vel(idim,iatom)+fcart_m(idim,iatom)*s)/scdot
     else
       vel_nexthalf(idim,iatom)=zero
     endif
   enddo
 enddo

!Computation of the next positions
 xcart_next(:,:)=xcart(:,:)+vel_nexthalf(:,:)*ab_mover%dtion

 if (zDEBUG)then
   write(std_out,*) 'a:',a
   write(std_out,*) 'b:',b
   write(std_out,*) 's:',s
   write(std_out,*) 'scdot:',scdot
 end if

!Convert back to xred (reduced coordinates)

 call xcart2xred(ab_mover%natom,rprimd,xcart_next,xred_next)

!write(std_out,*) 'isokinetic 06'
!##########################################################
!### 06. Update the history with the prediction

 xcart=xcart_next
 xred=xred_next

!increment the ihist
 hist%ihist = abihist_findIndex(hist,+1)

!Fill the history with the variables
!xred, acell, rprimd, vel
 call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 hist%vel(:,:,hist%ihist)=vel(:,:)
 hist%time(hist%ihist)=real(itime,kind=dp)*ab_mover%dtion

 if(zDEBUG)then
   write (std_out,*) 'fcart:'
   do kk=1,ab_mover%natom
     write (std_out,*) fcart(:,kk)
   end do
   write (std_out,*) 'vel:'
   do kk=1,ab_mover%natom
     write (std_out,*) vel(:,kk)
   end do
   write (std_out,*) 'strten:'
   write (std_out,*) strten(1:3),ch10,strten(4:6)
   write (std_out,*) 'etotal:'
   write (std_out,*) etotal
 end if

 if (.false.) write(std_out,*) ntime

end subroutine pred_isokinetic
!!***

end module m_pred_isokinetic
!!***

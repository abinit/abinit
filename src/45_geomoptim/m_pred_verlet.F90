!!****m* ABINIT/m_pred_verlet
!! NAME
!!  m_pred_verlet
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, JCC, SE)
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

module m_pred_verlet

 use defs_basis
 use m_abicore
 use m_abimover
 use m_abihist
 use m_xfpack

 use m_geometry,       only : mkrdim, xcart2xred, xred2xcart, fcart2fred, metric

 implicit none

 private
!!***

 public :: pred_verlet
!!***

contains
!!***

!!****f* ABINIT/pred_verlet
!! NAME
!! pred_verlet
!!
!! FUNCTION
!! Ionmov predictors (6 & 7) Verlet algorithm
!!
!! IONMOV 6:
!! Given a starting point xred that is a vector of length 3*natom
!! (reduced nuclei coordinates), a velocity vector (in cartesian
!! coordinates), and unit cell parameters (acell and rprimd -
!! without velocities in the present implementation), the Verlet
!! dynamics is performed, using the gradient of the energy
!! (atomic forces and stresses) as calculated by the routine scfcv.
!!
!! Some atoms can be kept fixed, while the propagation of unit cell
!! parameters is only performed if optcell/=0.
!! No more than ab_mover%ntime steps are performed.
!! The time step is governed by dtion, contained in ab_mover
!! (coming from dtset).
!! Returned quantities are xred, and eventually acell and rprimd (new ones!).
!!
!! IONMOV 7:
!! Block every atom for which the scalar product of velocity and
!! forces is negative, in order to reach the minimum.
!! The convergence requirement on the atomic forces, ab_mover%tolmxf,
!! allows an early exit.
!!
!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information needed by the preditor
!! itime  : Index of the present iteration
!! ntime  : Maximal number of iterations
!! ionmov : (6 or 7) Specific kind of VERLET
!! zDEBUG : if true print some debugging information
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces
!!                               acell, rprimd, stresses
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      fcart2fred,hist2var,metric,mkrdim,var2hist,wrtout,xcart2xred
!!      xfpack_f2vout,xfpack_vin2x,xfpack_x2vin,xred2xcart
!!
!! SOURCE

subroutine pred_verlet(ab_mover,hist,ionmov,itime,ntime,zDEBUG,iexit)

 implicit none

!Arguments ------------------------------------
!scalars
 type(abimover),intent(in) :: ab_mover
 type(abihist),intent(inout) :: hist
 integer,intent(in) :: itime
 integer,intent(in) :: ntime
 integer,intent(in) :: ionmov
 integer,intent(in) :: iexit
 logical,intent(in) :: zDEBUG

!Local variables-------------------------------
!scalars
 integer  :: ii,istopped,jj,kk,ndim,nstopped
 real(dp) :: amass_tot,etotal,diag,favg,ekin_corr,scprod,taylor,ucvol0
 real(dp),save :: ucvol,ucvol_next
 character(len=500) :: message
!arrays
 integer  :: stopped(ab_mover%natom)
 real(dp) :: acell0(3),fcart(3,ab_mover%natom)
 real(dp) :: fred_corrected(3,ab_mover%natom)
 real(dp) :: gprimd(3,3),gmet(3,3),rmet(3,3), strten(6)
 real(dp) :: xcart(3,ab_mover%natom),xcart_next(3,ab_mover%natom)
 real(dp) :: xred(3,ab_mover%natom),xred_next(3,ab_mover%natom)
 real(dp) :: vel(3,ab_mover%natom),vel_nexthalf(3,ab_mover%natom)
 real(dp),save :: acell(3),acell_next(3)
 real(dp),save :: rprimd(3,3),rprim(3,3),rprimd_next(3,3),rprim_next(3,3)
 real(dp),allocatable,save :: hessin(:,:)
 real(dp),allocatable,save :: vin(:),vin_prev(:),vin_next(:)
 real(dp),allocatable,save :: vout(:),vout_prev(:),vel_prevhalf(:,:)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 if(iexit/=0)then
   if (allocated(vin))           then
     ABI_DEALLOCATE(vin)
   end if
   if (allocated(vin_next))      then
     ABI_DEALLOCATE(vin_next)
   end if
   if (allocated(vout))          then
     ABI_DEALLOCATE(vout)
   end if
   if (allocated(vin_prev))      then
     ABI_DEALLOCATE(vin_prev)
   end if
   if (allocated(vout_prev))     then
     ABI_DEALLOCATE(vout_prev)
   end if
   if (allocated(hessin))        then
     ABI_DEALLOCATE(hessin)
   end if
   if (allocated(vel_prevhalf))  then
     ABI_DEALLOCATE(vel_prevhalf)
   end if
   return
 end if

!write(std_out,*) 'verlet 01'
!##########################################################
!### 01. Compute the dimension of vectors (ndim)

 ndim=3*ab_mover%natom
 if(ab_mover%optcell==1 .or.&
& ab_mover%optcell==4 .or.&
& ab_mover%optcell==5 .or.&
& ab_mover%optcell==6) ndim=ndim+1
 if(ab_mover%optcell==2 .or.&
& ab_mover%optcell==3) ndim=ndim+6
 if(ab_mover%optcell==7 .or.&
& ab_mover%optcell==8 .or.&
& ab_mover%optcell==9) ndim=ndim+3

!write(std_out,*) 'verlet 02'
!##########################################################
!### 02. Allocate the vectors vin, vout and hessian matrix

!Notice that vin, vout, etc could be allocated
!From a previous dataset with a different ndim
 if(itime==1)then
   if (allocated(vin))           then
     ABI_DEALLOCATE(vin)
   end if
   if (allocated(vin_next))      then
     ABI_DEALLOCATE(vin_next)
   end if
   if (allocated(vout))          then
     ABI_DEALLOCATE(vout)
   end if
   if (allocated(vin_prev))      then
     ABI_DEALLOCATE(vin_prev)
   end if
   if (allocated(vout_prev))     then
     ABI_DEALLOCATE(vout_prev)
   end if
   if (allocated(hessin))        then
     ABI_DEALLOCATE(hessin)
   end if
   if (allocated(vel_prevhalf))  then
     ABI_DEALLOCATE(vel_prevhalf)
   end if

   ABI_ALLOCATE(vin,(ndim))
   ABI_ALLOCATE(vin_next,(ndim))
   ABI_ALLOCATE(vout,(ndim))
   ABI_ALLOCATE(vin_prev,(ndim))
   ABI_ALLOCATE(vout_prev,(ndim))
   ABI_ALLOCATE(hessin,(ndim,ndim))
   ABI_ALLOCATE(vel_prevhalf,(3,ab_mover%natom))
 end if

!write(std_out,*) 'verlet 03'
!##########################################################
!### 03. Obtain the present values from the history

 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 do ii=1,3
   rprim(ii,1:3)=rprimd(ii,1:3)/acell(1:3)
 end do

 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)
 fcart(:,:)  =hist%fcart(:,:,hist%ihist)
 strten(:)  =hist%strten(:,hist%ihist)
 vel(:,:)   =hist%vel(:,:,hist%ihist)
 etotal     =hist%etot(hist%ihist)

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

 acell0=acell ; ucvol0=ucvol

!Get rid of mean force on whole unit cell, but only if no
!generalized constraints are in effect
 call fcart2fred(fcart,fred_corrected,rprimd,ab_mover%natom)
 if(ab_mover%nconeq==0)then
   amass_tot=sum(ab_mover%amass(:))
   do kk=1,3
     if (kk/=3.or.ab_mover%jellslab==0) then
       favg=sum(fred_corrected(kk,:))/dble(ab_mover%natom)
       fred_corrected(kk,:)=fred_corrected(kk,:)-favg*ab_mover%amass(:)/amass_tot
     end if
   end do
 end if

!write(std_out,*) 'verlet 04'
!##########################################################
!### 04. Fill the vectors vin and vout

!Initialize input vectors : first vin, then vout
 call xfpack_x2vin(acell, acell0, ab_mover%natom, ndim,&
& ab_mover%nsym, ab_mover%optcell, rprim, rprimd,&
& ab_mover%symrel, ucvol, ucvol0, vin, xred)
 call xfpack_f2vout(fred_corrected, ab_mover%natom, ndim,&
& ab_mover%optcell, ab_mover%strtarget, strten, ucvol,&
& vout)

!write(std_out,*) 'verlet 05'
!##########################################################
!### 05. Initialize or update the hessian matrix

!Here, set up the matrix of transformation between forces and
!acceleration. Masses must be included here.
!Beside this feature, one could define
!a preconditioner, in which case it should
!be the inverse hessian, like in Broyden. This explains the
!name chosen for this transformation matrix. This would allow
!to find easily the optimal geometry with ionmov=7.
!The default, now implemented, corresponds to the identity matrix
!in cartesian coordinates, which makes use of metric tensor gmet
!in reduced coordinates.

!Initialise the Hessian matrix using gmet
 if (itime==1)then
!  Initialize inverse hessian with identity matrix
!  in cartesian coordinates, which makes use of metric tensor gmet
!  in reduced coordinates.
   hessin(:,:)=zero
   do ii=1,ab_mover%natom
     do kk=1,3
       do jj=1,3
!        Warning : implemented in reduced coordinates
         if (ab_mover%iatfix(kk,ii)==0 .and.&
&         ab_mover%iatfix(jj,ii)==0 )then
           hessin(kk+3*(ii-1),jj+3*(ii-1))=gmet(kk,jj)/ab_mover%amass(ii)
         end if
       end do
     end do
   end do
   if(ab_mover%optcell/=0)then
!    These values might lead to too large changes in some cases ...
     diag=ab_mover%strprecon*30.0_dp/ucvol
     if(ab_mover%optcell==1) diag=diag/three
     do ii=3*ab_mover%natom+1,ndim
       hessin(ii,ii)=diag
     end do
   end if
 end if

!zDEBUG (vin,vout and hessin before prediction)
 if(zDEBUG)then
   write(std_out,*) 'Vectors vin and vout and inverse of Hessian (hessin) [before]'
   write(std_out,*) 'vin:'
   do ii=1,ndim,3
     if (ii+2<=ndim)then
       write(std_out,*) ii,vin(ii:ii+2)
     else
       write(std_out,*) ii,vin(ii:ndim)
     end if
   end do
   write(std_out,*) 'vout:'
   do ii=1,ndim,3
     if (ii+2<=ndim)then
       write(std_out,*) ii,vout(ii:ii+2)
     else
       write(std_out,*) ii,vout(ii:ndim)
     end if
   end do
   write(std_out,*) 'Inverse Hessian (hessin): ',ndim,'x',ndim
   do kk=1,ndim
     do jj=1,ndim,3
       if (jj+2<=ndim)then
         write(std_out,*) jj,hessin(jj:jj+2,kk)
       else
         write(std_out,*) jj,hessin(jj:ndim,kk)
       end if
     end do
   end do
 end if

!write(std_out,*) 'verlet 06'
!##########################################################
!### 06. Compute the next values

!%%% VERLET ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if(zDEBUG)then
   write(std_out,*) 'Shifted data ENTER'
   write(std_out,*) 'acell(:)',acell(:)
   write(std_out,*) 'rprimd(:,:)',rprimd(:,:)
   write(std_out,*) 'ucvol',ucvol
   write(std_out,*) 'vel_prevhalf(:,:)',vel_prevhalf(:,:)
   write(std_out,*) 'vin_prev(:)',vin_prev(:)
   write(std_out,*) 'vin(:)',vin(:)
   write(std_out,*) 'xcart(:,:)',xcart(:,:)
   write(std_out,*) 'xred(:,:)',xred(:,:)
 end if

!Compute next atomic coordinates and cell parameters, using
!Verlet algorithm

!1. First propagate the position, without acceleration
 if(itime/=1)then
   vin_next(:)=2*vin(:)-vin_prev(:)
   taylor=one
 else
!  Initialisation : no vin_prev is available, but the ionic velocity
!  is available, in cartesian coordinates
!  Uses the velocity
   xcart_next(:,:)=xcart(:,:)+ab_mover%dtion*vel(:,:)
!  Convert back to xred_next (reduced coordinates)
   call xcart2xred(ab_mover%natom,rprimd,xcart_next,xred_next)
!  Impose no change of acell, ucvol, rprim, and rprimd
   acell_next(:)=acell(:)
   ucvol_next=ucvol
   rprim_next(:,:)=rprim(:,:)
   rprimd_next(:,:)=rprimd(:,:)
   write(std_out,*) 'ucvol',ucvol
!  Store all these next values in vin_next
   call xfpack_x2vin(acell_next,acell0,ab_mover%natom,&
&   ndim,ab_mover%nsym,ab_mover%optcell,rprim_next,&
&   rprimd,ab_mover%symrel,ucvol_next,ucvol0,&
&   vin_next,xred_next)
   taylor=half
 end if

!2. Now, take into account the acceleration
 do ii=1,ndim
!  Note the minus sign: the forces are minus the gradients,
!  contained in vout.
   vin_next(:)=vin_next(:)-ab_mover%dtion**2*hessin(:,ii)*&
&   vout(ii)*taylor
 end do

!3. Implement fixing of atoms : put back old values for fixed
!components
 do kk=1,ab_mover%natom
   do jj=1,3
!    Warning : implemented in reduced coordinates
     if (ab_mover%iatfix(jj,kk) == 1) then
       vin_next(jj+(kk-1)*3)=vin(jj+(kk-1)*3)
     end if
   end do
 end do

!4. Now, compute the velocity at the next half-step
!Get xred_next, and eventually acell_next, ucvol_next, rprim_next and
!rprimd_next, from vin_next
 call xfpack_vin2x(acell_next,acell0,ab_mover%natom,ndim,&
& ab_mover%nsym,ab_mover%optcell,rprim_next,rprimd,&
& ab_mover%symrel,ucvol_next,ucvol0,vin_next,xred_next)
 if(ab_mover%optcell/=0)then
   call mkrdim(acell_next,rprim_next,rprimd_next)
   call metric(gmet,gprimd,-1,rmet,rprimd_next,ucvol_next)
 else
   rprimd_next(:,:)=rprimd(:,:)
 end if
!Convert input xred_next (reduced coordinates) to
!xcart_next (cartesian)
 call xred2xcart(ab_mover%natom,rprimd_next,xcart_next,xred_next)
!Compute the velocity at half of the new step
 vel_nexthalf(:,:)=(xcart_next(:,:)-xcart(:,:))/ab_mover%dtion

!5. If needed, compute the velocity at present position
 if(itime/=1)then
   vel(:,:)=(vel_nexthalf(:,:)+vel_prevhalf(:,:))*0.5_dp
 end if

!%%% VERLET ALGORITHM BLOCKING ATOMS %%%%%%%%%%%%%%%%%%%%%%

!Here, stop the atoms for which the scalar product of velocity
!and force is negative, and recompute the kinetic energy.
 if(ionmov==7)then
   stopped(:)=0
   do ii=1,ab_mover%natom
     scprod=fcart(1,ii)*vel(1,ii)+&
&     fcart(2,ii)*vel(2,ii)+&
&     fcart(3,ii)*vel(3,ii)
     if(scprod<0.0_dp .and. itime/=1)then
       stopped(ii)=1
       write(std_out,*) 'Stopped atom',ii
!      Shift the velocities of the previous half-step and current
!      half-step, so that the acceleration is correct but the
!      present velocity vanishes.
       vel_prevhalf(:,ii)=vel_prevhalf(:,ii)-vel(:,ii)
       vel_nexthalf(:,ii)=vel_nexthalf(:,ii)-vel(:,ii)
       vel(:,ii)=0.0_dp
       xcart_next(:,ii)=xcart(:,ii)+ab_mover%dtion*vel_nexthalf(:,ii)
     end if
   end do

   if(zDEBUG)then
     write (std_out,*) 'fcart:'
     do kk=1,ab_mover%natom
       write (std_out,*) fcart(:,kk)
     end do
     write (std_out,*) 'vel_prevhalf:'
     do kk=1,ab_mover%natom
       write (std_out,*) vel_prevhalf(:,kk)
     end do
     write (std_out,*) 'vel_nexthalf:'
     do kk=1,ab_mover%natom
       write (std_out,*) vel_nexthalf(:,kk)
     end do
     write (std_out,*) 'vel:'
     do kk=1,ab_mover%natom
       write (std_out,*) vel(:,kk)
     end do
     write (std_out,*) 'xcart_next:'
     do kk=1,ab_mover%natom
       write (std_out,*) xcart_next(:,kk)
     end do
   end if

!  Establish a list of stopped atoms
   nstopped=sum(stopped(:))

   if(nstopped/=0)then
     write(message,'(a)') ' List of stopped atoms (ionmov=7) :'
     call wrtout(ab_out,message,'COLL')
     istopped=1
     do ii=1,ab_mover%natom
       if(stopped(ii)==1)then
         stopped(istopped)=ii
         istopped=istopped+1
       end if
     end do
     do ii=1,nstopped,16
       write(message, '(16i4)' ) stopped(ii:min(ii+15,nstopped))
       call wrtout(ab_out,message,'COLL')
     end do
!    Now, compute the corrected kinetic energy
!    Generate xred_next from xcart_next
     call xcart2xred(ab_mover%natom,rprimd_next,xcart_next,xred_next)
!    Store xred_next, and eventual acell_next and rprim_next in vin
     call xfpack_x2vin(acell_next,acell0,&
&     ab_mover%natom,ndim,ab_mover%nsym,ab_mover%optcell,&
&     rprim_next,rprimd,&
&     ab_mover%symrel,ucvol_next,ucvol0,vin_next,xred_next)

     ekin_corr=0.0_dp
     do ii=1,ab_mover%natom
       do jj=1,3
!        Warning : the fixing of atomis is implemented in reduced
!        coordinates, so that this expression is wrong
         if (ab_mover%iatfix(jj,ii) == 0) then
           ekin_corr=ekin_corr+0.5_dp*ab_mover%amass(ii)*vel(jj,ii)**2
         end if
       end do
     end do
!    End of test nstopped/=0
   end if

!  End of test ionmov==7
 end if

!write(std_out,*) 'verlet 07'
!##########################################################
!### 07. Shift the data from next values to the present

!acell(:)=acell_next(:)
!rprim(:,:)=rprim_next(:,:)
!ucvol=ucvol_next
 vel_prevhalf(:,:)=vel_nexthalf(:,:)
 vin_prev(:)=vin(:)
 vin(:)=vin_next(:)
 xcart(:,:)=xcart_next(:,:)
 xred(:,:)=xred_next(:,:)

 write(std_out,*) 'Shifted data EXIT'
 write(std_out,*) 'acell(:)',acell(:)
 write(std_out,*) 'rprim(:,:)',rprim(:,:)
 write(std_out,*) 'rprimd(:,:)',rprimd(:,:)
 write(std_out,*) 'ucvol',ucvol
 write(std_out,*) 'vel_prevhalf(:,:)',vel_prevhalf(:,:)
 write(std_out,*) 'vin_prev(:)',vin_prev(:)
 write(std_out,*) 'vin(:)',vin(:)
 write(std_out,*) 'xred(:,:)',xred(:,:)


!write(std_out,*) 'verlet 08'
!##########################################################
!### 08. Update the history with the prediction

!Increase indexes
 hist%ihist=abihist_findIndex(hist,+1)

!Fill the history with the variables
!xred, acell, rprimd, vel
 call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 hist%vel(:,:,hist%ihist)=vel(:,:)
 hist%time(hist%ihist)=real(itime,kind=dp)*ab_mover%dtion

 if(zDEBUG)then
   write (std_out,*) 'vel:'
   do kk=1,ab_mover%natom
     write (std_out,*) vel(:,kk)
   end do
 end if

 if (.false.) write(std_out,*) ntime

end subroutine pred_verlet
!!***

end module m_pred_verlet
!!***

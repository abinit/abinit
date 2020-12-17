!!****m* ABINIT/m_pred_fire
!! NAME
!!  m_pred_fire
!!
!! FUNCTION
!! Ionmov predictors (15) FIRE algorithm
!! The fast inertial relaxation engine (FIRE) method for relaxation.
!! The method is described in  Erik Bitzek, Pekka Koskinen, Franz G"ahler, 
!! Michael Moseler, and Peter Gumbsch, Phys. Rev. Lett. 97, 170201 [[cite:Bitzek2006]]
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (hexu)
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

module m_pred_fire
 use defs_basis
 use m_abicore
 use m_abimover
 use m_abihist
 use m_xfpack
 use m_geometry,    only : mkrdim, fcart2fred, metric, xred2xcart
 use m_errors, only: unused_var

 implicit none

 private

 public :: pred_fire

contains

!!***

!!****f* m_pred_fire/pred_fire
!! NAME
!! pred_fire
!!
!! FUNCTION
!!
!! IONMOV 15:
!! Given a starting point xred that is a vector of length 3*natom
!! (reduced nuclei coordinates), a velocity vector (in cartesian
!! coordinates), and unit cell parameters (acell and rprimd -
!! without velocities in the present implementation), the Verlet
!! dynamics is performed, using the gradient of the energy
!! (atomic forces and stresses) as calculated by the routine scfcv.
!!
!! At each step, the dot product of velocity and force is calculated.
!! If the v.dot.f is positive for min_downhill consecutive steps, the
!! time step dtion is increased by dtinc until it reaches dtmax. If
!! v.dot.f is negative, dtion is mulitiplied by dtdec.
!!
!! Some atoms can be kept fixed, while the propagation of unit cell
!! parameters is only performed if optcell/=0.
!! No more than ab_mover%ntime steps are performed.
!! The time step is governed by dtion, contained in ab_mover
!! (coming from dtset).
!! Returned quantities are xred, and eventually acell and rprimd
!! (new ones!).
!!

!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information
!!                                needed by the preditor
!! itime  : Index of the present iteration
!! ntime  : Maximal number of iterations
!! ionmov : (15) FIRE. Not used in function. just to keep the same format as bfgs.
!! zDEBUG : if true print some debugging information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces
!!                               acell, rprimd, stresses
!!
!! PARENTS
!!      m_precpred_1geo
!!
!! CHILDREN
!!      fcart2fred,hist2var,metric,mkrdim,var2hist,xfpack_f2vout,xfpack_vin2x
!!      xfpack_x2vin,xred2xcart
!!
!! SOURCE
subroutine pred_fire(ab_mover, ab_xfh,forstr,hist,ionmov,itime,zDEBUG,iexit)

!Arguments ------------------------------------
!scalars
 type(abimover),intent(in) :: ab_mover
 type(ab_xfh_type),intent(inout) :: ab_xfh
 type(abiforstr),intent(in) :: forstr
 type(abihist),intent(inout) :: hist
 integer, intent(in) :: ionmov
 integer,intent(in) :: itime
 integer,intent(in) :: iexit
 logical,intent(in) :: zDEBUG


!Local variables-------------------------------
!scalars
integer  :: ihist,ihist_prev,ndim
integer, parameter :: min_downhill=4
integer  :: ii,jj,kk
real(dp),save :: ucvol0
real(dp) :: ucvol
real(dp) :: etotal,etotal_prev
real(dp) :: favg
! time step, damping factor initially dtion
real(dp),save :: dtratio, alpha
! dtinc: increment of dtratio
! dtdec: decrement of dtratio
! dtmax: maximum allowd value of dtratio
! alphadec: decrement of alpha
! alpha0: initial value of alpha
! mixold: if energy goes up, linear mix old and new coordinates. mixold
real(dp), parameter :: dtinc=1.1, dtdec=0.5, dtmax=10.0
real(dp), parameter :: alphadec=0.99, alpha0=0.2, mixold=0.3
! v.dot.f
real(dp) :: vf
! number of v.dot.f >0
integer, save :: ndownhill
! reset_lattice: whether to reset lattice if energy goes up.
logical, parameter :: reset_lattice = .true.

!arrays
real(dp) :: gprimd(3,3)
real(dp) :: gmet(3,3)
real(dp) :: rmet(3,3)
real(dp) :: fcart(3, ab_mover%natom)
real(dp) :: acell(3),strten(6), acell0(3)
real(dp) :: rprim(3,3),rprimd(3,3), rprimd0(3,3)
real(dp) :: xred(3,ab_mover%natom),xcart(3,ab_mover%natom)
! velocity are saved
real(dp) :: vel(3,ab_mover%natom)
real(dp) :: residual(3,ab_mover%natom),residual_corrected(3,ab_mover%natom)
real(dp),allocatable, save :: vin(:), vout(:)
real(dp), allocatable, save:: vin_prev(:)
! velocity but correspoing to vin&vout, for ion&cell relaxation
real(dp),allocatable,save :: vel_ioncell(:)

 ABI_UNUSED((/ionmov, ab_xfh%mxfh/))

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 if(iexit/=0)then
   if (allocated(vin))           then
     ABI_DEALLOCATE(vin)
   end if
   if (allocated(vout))          then
     ABI_DEALLOCATE(vout)
   end if
   if (allocated(vin_prev))           then
     ABI_DEALLOCATE(vin_prev)
   end if
   if (allocated(vel_ioncell))          then
     ABI_DEALLOCATE(vel_ioncell)
   end if
   return
 end if

 !write(std_out,*) 'FIRE 01'
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

!write(std_out,*) 'FIRE: ndim=', ndim
! write(std_out,*) 'FIRE 02'
!##########################################################
!### 02. Allocate the vectors vin

!Notice that vin, vout, etc could be allocated
!From a previous dataset with a different ndim
 if(itime==1)then
   if (allocated(vin))           then
     ABI_DEALLOCATE(vin)
   end if
   if (allocated(vout))          then
     ABI_DEALLOCATE(vout)
   end if
   if (allocated(vin_prev))           then
     ABI_DEALLOCATE(vin_prev)
   end if
   if (allocated(vel_ioncell))          then
     ABI_DEALLOCATE(vel_ioncell)
   end if

   ABI_ALLOCATE(vin,(ndim))
   ABI_ALLOCATE(vout,(ndim))
   ABI_ALLOCATE(vin_prev,(ndim))
   ABI_ALLOCATE(vel_ioncell,(ndim))
   vel_ioncell(:)=0.0
 end if

 !write(std_out,*) 'FIRE 03'
!##########################################################
!### 03. Obtain the present values from the history

 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 do ii=1,3
   rprim(ii,1:3)=rprimd(ii,1:3)/acell(1:3)
 end do

 ihist = abihist_findIndex(hist, 0)
 ihist_prev  = abihist_findIndex(hist,-1)

 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)
 fcart(:,:)  =hist%fcart(:,:,hist%ihist)
 strten(:)  =hist%strten(:,hist%ihist)
 etotal=hist%etot(hist%ihist)

 if(itime==1) then
     etotal_prev=0.0
 else
     etotal_prev=hist%etot(ihist_prev)
 endif

!Fill the residual with forces (No preconditioning)
!Or the preconditioned forces
 if (ab_mover%goprecon==0)then
   call fcart2fred(hist%fcart(:,:,hist%ihist),residual,rprimd,ab_mover%natom)
 else
   residual(:,:)=forstr%fred(:,:)
 end if

!Save initial values
 if (itime==1)then
   acell0(:)=acell(:)
   rprimd0(:,:)=rprimd(:,:)
   ucvol0=ucvol
 end if


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

 residual_corrected(:,:)=residual(:,:)
 if(ab_mover%nconeq==0)then
   do kk=1,3
     if (kk/=3.or.ab_mover%jellslab==0) then
       favg=sum(residual_corrected(kk,:))/dble(ab_mover%natom)
       residual_corrected(kk,:)=residual_corrected(kk,:)-favg
     end if
   end do
 end if

 !write(std_out,*) 'FIRE 04'
!##########################################################
!### 04. Fill the vectors vin and vout

!Initialize input vectors : first vin, then vout
! transfer xred, acell, and rprim to vin
call xfpack_x2vin(acell, acell0, ab_mover%natom, ndim,&
& ab_mover%nsym, ab_mover%optcell, rprim, rprimd0,&
& ab_mover%symrel, ucvol, ucvol0, vin, xred)
!end if

!transfer fred and strten to vout. 
!Note: fred is not f in reduced co.
!but dE/dx

 call xfpack_f2vout(residual_corrected, ab_mover%natom, ndim,&
& ab_mover%optcell, ab_mover%strtarget, strten, ucvol,&
& vout)
! Now vout -> -dE/dx
vout(:) = -1.0*vout(:)

 !write(std_out,*) 'FIRE 05'
!##########################################################
!### 05. iniialize FIRE
if ( itime==1 ) then
   ndownhill=0
   alpha=alpha0
   if (ab_mover%dtion>0)then
     dtratio = 1.0
   end if
end if

 !write(std_out,*) 'FIRE 06'
!##########################################################
!### 06. update timestep
! Note that vin & vout are in reduced coordinates.
vf=sum(vel_ioncell*vout)
if ( vf >= 0.0_dp .and. (etotal- etotal_prev <0.0_dp) ) then
!if ( vf >= 0.0_dp ) then
    ndownhill=ndownhill+1
    ! mix v with the v projected on force vector.
    vel_ioncell(:)=(1.0-alpha)*vel_ioncell(:) + alpha* vout *  &
&               sqrt(sum(vel_ioncell*vel_ioncell)/sum(vout*vout))
    if ( ndownhill>min_downhill ) then
        dtratio = min(dtratio * dtinc, dtmax)
        alpha = alpha * alphadec
    end if
else
    ! reset downhill counter, velocity, alpha. decrease dtratio.
    ndownhill=0
    vel_ioncell(:)=0.0
    alpha=alpha0
    dtratio = dtratio*dtdec
endif

 !write(std_out,*) 'FIRE 07'
!##########################################################
!### 07. MD step. update vel_ioncell

! Here mass is not used: all masses=1
!write(std_out,*) 'FIRE vin: ', vin
! update v
vel_ioncell = vel_ioncell + dtratio*ab_mover%dtion* vout
!write(std_out,*) 'FIRE vel: ', vel_ioncell
!write(std_out,*) 'FIRE delta x',dtratio*ab_mover%dtion* vel_ioncell
! update x
vin = vin + dtratio*ab_mover%dtion* vel_ioncell
!write(std_out,*) 'FIRE vin: ', vin
   
!   write(std_out,*) 'FIRE vout: ', vout
!   write(std_out,*) 'FIRE vf: ', vf
!   write(std_out,*) 'FIRE etotal: ', etotal
!   write(std_out,*) 'FIRE etotal_prev: ', etotal_prev
!   write(std_out,*) 'FIRE deltaE: ',etotal-etotal_prev
   write(std_out,*) 'FIRE ndownhill: ', ndownhill
!   write(std_out,*) 'FIRE dtratio: ', dtratio
!   write(std_out,*) 'FIRE dtion: ', ab_mover%dtion


!Implement fixing of atoms : put back old values for fixed
!components
 do kk=1,ab_mover%natom
   do jj=1,3
!    Warning : implemented in reduced coordinates
     if ( ab_mover%iatfix(jj,kk)==1) then
       vin(jj+(kk-1)*3)=vin_prev(jj+(kk-1)*3)
     end if
   end do
 end do

! reset_lattice to last step by a ratio if energy is increased.
! disabled for debugging
if ( etotal - etotal_prev >0.0 ) then
    vin= vin*(1-mixold)+vin_prev*mixold
end if

! only set vin to vin_prev when energy decreased, so it's 
! possible to go back.
! if (etotal - etotal_prev <0.0 ) then
 vin_prev(:)=vin(:)
! endif



!##########################################################
!### 08. update hist.

!Increase indexes
 hist%ihist = abihist_findIndex(hist,+1)
!Transfer vin  to xred, acell and rprim
 call xfpack_vin2x(acell, acell0, ab_mover%natom, ndim,&
& ab_mover%nsym, ab_mover%optcell, rprim, rprimd0,&
& ab_mover%symrel, ucvol, ucvol0,&
& vin, xred)


 if(ab_mover%optcell/=0)then
   call mkrdim(acell,rprim,rprimd)
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 end if

!Fill the history with the variables
!xcart, xred, acell, rprimd
 call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 ihist_prev = abihist_findIndex(hist,-1)
 hist%vel(:,:,hist%ihist)=hist%vel(:,:,ihist_prev)

 if(zDEBUG)then
   write (std_out,*) 'residual:'
   do kk=1,ab_mover%natom
     write (std_out,*) residual(:,kk)
   end do
   write (std_out,*) 'strten:'
   write (std_out,*) strten(1:3),ch10,strten(4:6)
   write (std_out,*) 'etotal:'
   write (std_out,*) etotal
 end if


end subroutine pred_fire
!!***

end module m_pred_fire


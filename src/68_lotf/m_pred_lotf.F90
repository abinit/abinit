!!****m* ABINIT/m_pred_lotf
!! NAME
!! m_pred_lotf
!!
!! FUNCTION
!! Contains the predictor for LOTF (ionmov==23)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, JCC, SE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module m_pred_lotf

 use m_abicore
 use defs_basis
 use m_abimover
 use m_abihist
 use lotfpath

 use m_numeric_tools, only : uniformrandom
 use m_geometry,      only : xcart2xred

 implicit none

 private

 public :: pred_lotf

CONTAINS !===========================================================
 !!***

 !!****f* ABINIT/m_pred_lotf/pred_lotf
 !! NAME
 !! pred_lotf
 !!
 !! FUNCTION
 !! Ionmov predictors (12) Lotf ensemble molecular dynamics
 !!
 !! IONMOV 23:
 !! Lotf ensemble molecular dynamics.
 !!
 !! COPYRIGHT
 !! Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, JCC, SE)
 !! This file is distributed under the terms of the
 !! GNU General Public License, see ~abinit/COPYING
 !! or http://www.gnu.org/copyleft/gpl.txt .
 !! For the initials of contributors,
 !! see ~abinit/doc/developers/contributors.txt .
 !!
 !! INPUTS
 !! ab_mover <type(abimover)> : Datatype with all the information
 !!                                needed by the preditor
 !! itime  : Index of the present iteration
 !! icycle : Index of the cycle in the iteration
 !! zDEBUG : if true print some debugging information
 !!
 !! OUTPUT
 !!
 !! SIDE EFFECTS
 !! hist <type(abihist)> : History of positions,forces
 !!                               acell, rprimd, stresses
 !!
 !! NOTES
 !!
 !! PARENTS
!!      mover
!!
 !! CHILDREN
!!      extrapolation_loop,fitclus,hist2var,init_lotf,intparms
!!      lotf_interpolation,var2hist,vel_rescale,vel_to_gauss,wrtout,xcart2xred
!!      xred2xcart
!!
 !! SOURCE
 subroutine pred_lotf(ab_mover,hist,itime,icycle,zDEBUG,iexit)

 use m_geometry,       only : xred2xcart
  implicit none

  !Arguments ------------------------
  type(abimover),intent(in)       :: ab_mover
  type(abihist),intent(inout) :: hist
  integer,intent(in) :: itime
  integer,intent(in) :: icycle
  integer,intent(in) :: iexit
  logical,intent(in) :: zDEBUG

  !Local variables-------------------------------
  !scalars
  integer  :: kk,iatom,idim,idum=5
  real(dp) :: v2gauss
  real(dp),parameter :: v2tol=tol8
  real(dp) :: etotal
  character(len=5000) :: message
  logical :: do_extrap,do_interp,do_first
  !arrays
  real(dp),allocatable,save :: fcart_m(:,:),vel_nexthalf(:,:)
  real(dp),allocatable,save :: xcart_old(:,:),vel_old(:,:)

  real(dp) :: acell(3),rprimd(3,3),fcart(3,ab_mover%natom)
  real(dp) :: xcart(3,ab_mover%natom),xcart_next(3,ab_mover%natom)
  real(dp) :: xred(3,ab_mover%natom),xred_next(3,ab_mover%natom)
  real(dp) :: vel(3,ab_mover%natom)
  real(dp) :: strten(6)

  !***************************************************************************
  !Beginning of executable session
  !***************************************************************************

  if(iexit/=0)then
    if (allocated(fcart_m))       then
      ABI_DEALLOCATE(fcart_m)
    end if
    if (allocated(vel_nexthalf))  then
      ABI_DEALLOCATE(vel_nexthalf)
    end if
    if (allocated(xcart_old))       then
      ABI_DEALLOCATE(xcart_old)
    end if
    if (allocated(vel_old))       then
      ABI_DEALLOCATE(vel_old)
    end if
    return
  end if

  !write(std_out,*) 'lotf 01'
  !##########################################################
  !### 01. Debugging and Verbose

  if(zDEBUG)then
    write(std_out,'(a,3a,40a,37a)') ch10,('-',kk=1,3),&
&   'Debugging and Verbose for pred_lotf',('-',kk=1,37)
    write(std_out,*) 'ionmov: ',12
    write(std_out,*) 'itime:  ',itime
  end if

  !write(std_out,*) 'lotf 02'
  !##########################################################
  !### 02. Allocate the vectors vin, vout and hessian matrix
  !###     These arrays could be allocated from a previus
  !###     dataset that exit before itime==ntime


  if (.not.allocated(fcart_m))       then
    ABI_ALLOCATE(fcart_m,(3,ab_mover%natom))
  end if
  if (.not.allocated(vel_nexthalf))  then
    ABI_ALLOCATE(vel_nexthalf,(3,ab_mover%natom))
  end if
  if (.not.allocated(vel_old))       then
    ABI_ALLOCATE(vel_old,(3,ab_mover%natom))
  end if
  if (.not.allocated(xcart_old))  then
    ABI_ALLOCATE(xcart_old,(3,ab_mover%natom))
  end if


  !write(std_out,*) 'lotf 03'
  !##########################################################
  !### 03. Set what Pred_lotf has to do

  !--First step only the first time pred_lotf is called
  do_first = (itime==1 .and. icycle==1)

  !--Extrapolation is computed only is itime is a multiple of lotfvar%nitex
  do_extrap = lotf_extrapolation(itime)

  !--Interpolation is done at icycle==2 when extrapolation is active,
  !--else at the first cycle
  if(do_extrap) then
    do_interp = (icycle==2)
  else
    do_interp = (icycle==1)
  endif

  !write(std_out,*) 'lotf 04'
  !##########################################################
  !### 04. Obtain the present values from the history

  call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)

  fcart(:,:)=hist%fcart(:,:,hist%ihist)
  strten(:) =hist%strten(:,hist%ihist)
  vel(:,:)  =hist%vel(:,:,hist%ihist)
  etotal    =hist%etot(hist%ihist)

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

  !write(std_out,*) 'lotf 05'
  !##########################################################
  !### 05. First half-time (First cycle the loop is double)

  first_step: if(do_first) then
    write(message, '(a)' ) ' --- LOTF INITIALIZATION---'
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')

    !--convert from xred to xcart
    call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

    !PRINT *,'HH-first-1',itime,icycle,sum(abs(xcart(1,:))),sum(abs(fcart(1,:))),sum(abs(vel(1,:)))

    !--call the LOTF initialization
    call init_lotf(itime,ab_mover%natom,acell,rprimd,xcart)

    !--Application of Gauss' principle of least constraint according to 
    ! Fei Zhang's algorithm (J. Chem. Phys. 106, 1997, p.6102 [[cite:Zhang1997]])
    !--v2gauss is twice the kinetic energy
    call vel_to_gauss(vel,ab_mover%amass,v2gauss)

    !--If there is no kinetic energy
    if(v2gauss<=v2tol) then
      !--Maxwell-Boltzman distribution
      v2gauss=zero
      do iatom=1,ab_mover%natom
        do idim=1,3
          vel(idim,iatom)=sqrt(kb_HaK*ab_mover%mdtemp(1)/ab_mover%amass(iatom))*cos(two_pi*uniformrandom(idum))
          vel(idim,iatom)=vel(idim,iatom)*sqrt(-2._dp*log(uniformrandom(idum)))
        end do
      end do
    endif

    !--Init alpha_end parameters
    call fitclus(.true.,fcart(:,:),xcart,alpha_end,1,ab_mover%natom)

    !--Forces/mass for the first step
    forall(iatom = 1:ab_mover%natom) fcart_m(:,iatom) = fcart(:,iatom)/ab_mover%amass(iatom)

    !print *,'b-end 1step',itime,sum(sum(alpha,dim=2)),sum(sum(alpha_in,dim=2)),sum(sum(alpha_end,dim=2))
    !PRINT *,'HH-first-2',itime,icycle,sum(abs(xcart(1,:))),sum(abs(fcart(1,:))),sum(abs(vel(1,:)))

  end if first_step


  !--Do extrapolation:
  !  when extrapolation is active for a fixed itime, mover procedure
  !  will call pred_lotf twice,
  !  1-a) The final values of alpha are set as initial value
  !  1-b) The actual coordinates and velocities are stored in a local saved array
  !  1-c) I the first the coordinated of the atoms at the step itime+lotfvar%nitex
  !        will be extrapolated and alpha_in updated
  !  2-a) The new final value of alpha_end (step itime+lotfvar%nitex) are computed
  !  2-b) Restore velocities and positions
  extrapolation_lotf: if(do_extrap) then
    select case(icycle)
    case(1)
      !--
      write(message, '(a)' )' LOTF : Extrapolation'
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')

      !--Intial alpha values to do extrpolation correspond the last one
      alpha_in = alpha_end

      !--store the value of xcart
      xcart_old = xcart
      vel_old = vel

    !PRINT *,'HH-extra-1',itime,icycle,sum(abs(xcart(1,:))),sum(abs(fcart(1,:))),sum(abs(vel(1,:)))

      !--Compute new bond, and positions by extrpolation (xcartfit) at
      !  the step=itime+lotfvar%nitex
      call extrapolation_loop(itime,ab_mover%mdtemp(1),ab_mover%dtion,ab_mover%amass,xcart,&
&                     vel,xcart_next,vel_nexthalf)

      !print *,'b-end extra',itime,sum(sum(alpha,dim=2)),sum(sum(alpha_in,dim=2)),sum(sum(alpha_end,dim=2))

      !--Convert back to xred (reduced coordinates)
      call xcart2xred(ab_mover%natom,rprimd,xcart_next,xred_next)

    !PRINT *,'HH-extra-2',itime,icycle,sum(abs(xcart(1,:))),sum(abs(fcart(1,:))),sum(abs(vel(1,:)))

    case(2)
    !PRINT *,'HH-extra-3',itime,icycle,sum(abs(xcart(1,:))),sum(abs(fcart(1,:))),sum(abs(vel(1,:)))

      call fitclus(.true.,fcart(:,:),xcart,alpha_end,1,ab_mover%natom)
      !print *,'b-end fitclus',itime,sum(sum(alpha,dim=2)),sum(sum(alpha_in,dim=2)),sum(sum(alpha_end,dim=2))

      !--restore the value of xcart
      xcart = xcart_old
      vel = vel_old
    !PRINT *,'HH-extra-4',itime,icycle,sum(abs(xcart(1,:))),sum(abs(fcart(1,:))),sum(abs(vel(1,:)))

    end select
  endif extrapolation_lotf


  !--Interpolation(pass here for extrapolation when icycle==2 else when icycle==1)
  if(do_interp)then

    !PRINT *,'HH-norma-1',itime,icycle,sum(abs(xcart(1,:))),sum(abs(fcart(1,:))),sum(abs(vel(1,:)))

    !--Compute rescaled vel (center of mass speed and gaussian)
    call vel_rescale(ab_mover%mdtemp(1),vel,ab_mover%amass,v2gauss)

    !PRINT *,'HH-norma-2',itime,icycle,sum(abs(xcart(1,:))),sum(abs(fcart(1,:))),sum(abs(vel(1,:)))

    !--VERLETVEL (interpolation) using the interpolated value of alpha
    !  in the interval [alpha_in,alpha_end]
    call lotf_interpolation(itime,ab_mover%dtion,v2gauss,ab_mover%amass,xcart,vel,fcart_m,xcart_next,vel_nexthalf)

    !PRINT *,'HH-norma-3',itime,icycle,sum(abs(xcart(1,:))),sum(abs(fcart(1,:))),sum(abs(vel(1,:)))

    !print *,'b-end interpo',itime,sum(sum(alpha,dim=2)),sum(sum(alpha_in,dim=2)),sum(sum(alpha_end,dim=2))

    !--Convert back to xred (reduced coordinates)
    call xcart2xred(ab_mover%natom,rprimd,xcart_next,xred_next)

    !--Interpolate alpha_in,alpha_end to obtain: alpha for the next step
    call intparms(itime)
    !print *,'b-end intparm',itime,sum(sum(alpha,dim=2)),sum(sum(alpha_in,dim=2)),sum(sum(alpha_end,dim=2))
  endif

  !write(std_out,*) 'lotf 06'
  !##########################################################
  !### 06. Update the history with the prediction

  xcart = xcart_next
  xred = xred_next

  !PRINT *,'HH-final-1',itime,icycle,sum(abs(xcart(1,:))),sum(abs(fcart(1,:))),sum(abs(vel(1,:)))

  !print *,"XCART_final",sum(abs(xcart(1,:))),sum(abs(xcart(2,:))),sum(abs(xcart(3,:)))

  write(message, '(a,3d20.9)' )&
&    ' --- XCART_final',sum(abs(xcart(1,:))),sum(abs(xcart(2,:))),sum(abs(xcart(3,:)))
  call wrtout(std_out,message,'PERS')

  !Increase indexes
  hist%ihist = abihist_findIndex(hist,+1)

  !Fill the history with the variables
  !xcart, xred, acell, rprimd
  call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
  hist%vel(:,:,hist%ihist) = vel(:,:)
  hist%time(hist%ihist)=real(itime,kind=dp)*ab_mover%dtion

  if(zDEBUG)then
    write (std_out,*) 'xcart:'
    do kk=1,ab_mover%natom
      write (std_out,*) xcart(:,kk)
    end do
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

 end subroutine pred_lotf
 !!***
end module m_pred_lotf
!!***

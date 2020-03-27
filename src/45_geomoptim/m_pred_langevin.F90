!!****m* ABINIT/m_pred_langevin
!! NAME
!!  m_pred_langevin
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

module m_pred_langevin

 use defs_basis
 use m_abicore
 use m_abimover
 use m_abihist

 use m_numeric_tools,  only : uniformrandom
 use m_geometry,    only : xcart2xred, xred2xcart, metric
 use m_results_gs , only : results_gs_type

 implicit none

 private
!!***

 public :: pred_langevin
!!***

contains
!!***

!!****f* ABINIT/pred_langevin
!! NAME
!! pred_langevin
!!
!! FUNCTION
!! Ionmov predictors (9) Langevin dynamics algorithm
!!
!! IONMOV 9:
!! Uses a Langevin dynamics algorithm :
!! see J. Chelikowsky, J. Phys. D : Appl Phys. 33(2000)R33 [[cite:Chelikowsky2000]]
!!
!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information
!!                                needed by the preditor
!! itime  : Index of the present iteration
!! ntime  : Maximal number of iterations
!! icycle : Index of the present cycle
!! ncycle : Maximal number of cycles
!! zDEBUG : if true print some debugging information
!!
!! OUTPUT
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

subroutine pred_langevin(ab_mover,hist,icycle,itime,ncycle,ntime,zDEBUG,iexit,skipcycle)

 implicit none

!Arguments ------------------------------------
!scalars
 type(abimover),intent(in)       :: ab_mover
 type(abihist),intent(inout) :: hist
 integer,intent(in)    :: itime
 integer,intent(in)    :: ntime
 integer,intent(in)    :: iexit
 integer,intent(in)    :: icycle
 integer,intent(inout) :: ncycle
 logical,intent(in)    :: zDEBUG
 logical,intent(out)   :: skipcycle

!Local variables-------------------------------
!scalars
 integer  :: ii,kk,iatom,idim,iatom1,iatom2,itypat,idum=5,ihist_prev,mcfac
 real(dp) :: ucvol,ucvol_next
 real(dp),parameter :: v2tol=tol8
 real(dp) :: etotal,rescale_vel,ran_num1,ran_num2
 real(dp) :: ktemp,dist,distx,disty,distz,maxp1,maxp2,v2nose
 real(dp) :: sig_gauss,delxi
 logical  :: jump_end_of_cycle=.FALSE.
 character(len=5000) :: message
!arrays
 integer,allocatable :: imax_perm(:)
 real(dp),allocatable,save :: max_perm(:),pot_perm(:)
 real(dp),allocatable,save :: ran_force(:,:),lang_force(:,:)
 real(dp),allocatable,save :: fcart_mold(:,:),fcart_m(:,:)

 real(dp) :: acell(3),acell_next(3)
 real(dp) :: rprim(3,3),rprimd(3,3),rprimd_next(3,3),rprim_next(3,3)
 real(dp) :: gprimd(3,3),gmet(3,3),rmet(3,3)
 real(dp) :: fcart(3,ab_mover%natom)
!real(dp) :: fred_corrected(3,ab_mover%natom)
 real(dp) :: xcart(3,ab_mover%natom),xcart_next(3,ab_mover%natom)
 real(dp) :: xred(3,ab_mover%natom),xred_next(3,ab_mover%natom)
 real(dp) :: vel(3,ab_mover%natom)
 real(dp) :: strten(6)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 if(iexit/=0)then
   if (allocated(pot_perm))    then
     ABI_DEALLOCATE(pot_perm)
   end if
   if (allocated(max_perm))    then
     ABI_DEALLOCATE(max_perm)
   end if
   if (allocated(imax_perm))   then
     ABI_DEALLOCATE(imax_perm)
   end if
   if (allocated(ran_force))   then
     ABI_DEALLOCATE(ran_force)
   end if
   if (allocated(lang_force))  then
     ABI_DEALLOCATE(lang_force)
   end if
   if (allocated(fcart_mold))  then
     ABI_DEALLOCATE(fcart_mold)
   end if
   if (allocated(fcart_m))     then
     ABI_DEALLOCATE(fcart_m)
   end if
   return
 end if

 jump_end_of_cycle=.FALSE.

!write(std_out,*) 'langevin 02',jump_end_of_cycle
!##########################################################
!### 02. Allocate the arrays
!###     These arrays could be allocated from a previus
!###     dataset that exit before itime==ntime

 if(itime==1)then
   if (allocated(pot_perm))    then
     ABI_DEALLOCATE(pot_perm)
   end if
   if (allocated(max_perm))    then
     ABI_DEALLOCATE(max_perm)
   end if
   if (allocated(imax_perm))   then
     ABI_DEALLOCATE(imax_perm)
   end if
   if (allocated(ran_force))   then
     ABI_DEALLOCATE(ran_force)
   end if
   if (allocated(lang_force))  then
     ABI_DEALLOCATE(lang_force)
   end if
   if (allocated(fcart_mold))  then
     ABI_DEALLOCATE(fcart_mold)
   end if
   if (allocated(fcart_m))     then
     ABI_DEALLOCATE(fcart_m)
   end if
 end if

 if (.not.allocated(pot_perm))    then
   ABI_ALLOCATE(pot_perm,(ab_mover%natom))
 end if
 if (.not.allocated(max_perm))    then
   ABI_ALLOCATE(max_perm,(ab_mover%ntypat))
 end if
 if (.not.allocated(imax_perm))   then
   ABI_ALLOCATE(imax_perm,(ab_mover%ntypat))
 end if
 if (.not.allocated(ran_force))   then
   ABI_ALLOCATE(ran_force,(3,ab_mover%natom))
 end if
 if (.not.allocated(lang_force))  then
   ABI_ALLOCATE(lang_force,(3,ab_mover%natom))
 end if
 if (.not.allocated(fcart_mold))  then
   ABI_ALLOCATE(fcart_mold,(3,ab_mover%natom))
 end if
 if (.not.allocated(fcart_m))     then
   ABI_ALLOCATE(fcart_m,(3,ab_mover%natom))
 end if

!write(std_out,*) 'langevin 03',jump_end_of_cycle
!##########################################################
!### 03. Obtain the present values from the history

 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)

 fcart(:,:)=hist%fcart(:,:,hist%ihist)
 strten(:) =hist%strten(:,hist%ihist)
 vel(:,:)  =hist%vel(:,:,hist%ihist)
 etotal    =hist%etot(hist%ihist)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 do ii=1,3
   rprim(ii,1:3)=rprimd(ii,1:3)/acell(1:3)
 end do
 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

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

!write(std_out,*) 'langevin 04',jump_end_of_cycle
!##########################################################
!### 04. Compute the next values (Only for the first cycle)

 if (icycle==1) then

!  The temperature is linear between initial and final values
!  It is here converted from Kelvin to Hartree (kb_HaK)
   ktemp=(ab_mover%mdtemp(1)+((ab_mover%mdtemp(2)-ab_mover%mdtemp(1))/dble(ntime-1))*(itime-1))*kb_HaK
!  write(std_out,*) 'KTEMP=',ktemp
!  write(std_out,*) 'MDITEMP=',ab_mover%mdtemp(1)
!  write(std_out,*) 'MDFTEMP=',ab_mover%mdtemp(2)
!  write(std_out,*) 'DELAYPERM=',ab_mover%delayperm


!  %%% LANGEVIN DYNAMICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  1. Next values are the present ones
   acell_next(:)=acell(:)
   ucvol_next=ucvol
   rprim_next(:,:)=rprim(:,:)
   rprimd_next(:,:)=rprimd(:,:)

   if (zDEBUG) then
     write (std_out,*) '1. Next values are the present ones'
     write(std_out,*) 'RPRIMD'
     do kk=1,3
       write(std_out,*) rprimd(:,kk)
     end do
     write(std_out,*) 'RPRIM'
     do kk=1,3
       write(std_out,*) rprim(:,kk)
     end do
     write(std_out,*) 'ACELL'
     write(std_out,*) acell(:)
   end if

   if(itime==1)then

!    2.   Compute twice the kinetic energy of the system, called v2nose
     v2nose=0.0_dp
     do iatom=1,ab_mover%natom
       do idim=1,3
         v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*ab_mover%amass(iatom)
       end do
     end do
     if (zDEBUG) then
       write (std_out,*) '2.   Compute twice the kinetic energy of the system, called v2nose'
       write (std_out,*) 'V2NOSE=',v2nose
     end if

!    3.   If there is no kinetic energy, use random numbers
     if (v2nose<=v2tol) then
       v2nose=0.0_dp
       do iatom=1,ab_mover%natom
         do idim=1,3
!          uniformrandom returns a uniform random deviate between 0.0 and 1.0
!          if it were always 0 or 1, then the following expression
!          would give the requested temperature
           vel(idim,iatom)=(1.0_dp-2.0_dp*uniformrandom(idum))*sqrt(ktemp/ab_mover%amass(iatom))
!          Recompute v2nose
           v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*ab_mover%amass(iatom)
         end do
       end do
     end if

     if (zDEBUG) then
       write (std_out,*) '3.   If there is no kinetic energy, use random numbers'
       write (std_out,*) 'VEL'
       do kk=1,ab_mover%natom
         write (std_out,*) vel(:,kk)
       end do
       write (std_out,*) 'V2NOSE=',v2nose
     end if


!    Now, rescale the velocities to give the proper temperature
     rescale_vel=sqrt(3.0_dp*ab_mover%natom*(ab_mover%mdtemp(1))*kb_HaK/v2nose)
     vel(:,:)=vel(:,:)*rescale_vel
!    Recompute v2nose with the rescaled velocities
     v2nose=0.0_dp
     do iatom=1,ab_mover%natom
       do idim=1,3
         v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*ab_mover%amass(iatom)
       end do
     end do
     write(message, '(a)' )&
&     ' Rescaling or initializing velocities to initial temperature'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     write(message, '(a,D12.5,a,D12.5)' )&
&     ' ---  Scaling factor : ',rescale_vel,' Asked T (K) ',ab_mover%mdtemp(1)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     write(message, '(a,D12.5)' )&
&     ' ---  Effective temperature',v2nose/3.0_dp/(kb_HaK*ab_mover%natom)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
!    end if itime==0
   end if

!  This section is devoted to the optional atom permutation (JYR 001114)
!  Two input variables are needed
!  ab_mover%delayperm : is the interval (in time steps) at which
!  atoms are tentatively permuted
!  default value could be 0
!  ab_mover%signperm  : is the type of bias for the permutation
!  +1  to favor alternation of species
!  -1  to favor segregation

!  Force no permutation at initial step
   if (itime/=1 .and. ab_mover%delayperm/=0 .and. ab_mover%ntypat>2) then
     if (mod(itime-1,ab_mover%delayperm)==0) then
!      Try commutation of atoms.
       write(message, '(a)')' Attempt of commutation '
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
!      Compute a 'permutation potential'
       do iatom=1,ab_mover%natom
         pot_perm(iatom)=0.0_dp
         do iatom1=1,ab_mover%natom
           if (iatom1.ne.iatom) then
             distx=xcart(1,iatom)-xcart(1,iatom1)
             distx=distx-acell(1)*nint(distx/acell(1))
             disty=xcart(2,iatom)-xcart(2,iatom1)
             disty=disty-acell(2)*nint(disty/acell(2))
             distz=xcart(3,iatom)-xcart(3,iatom1)
             distz=distz-acell(3)*nint(distz/acell(3))
!            Here we count each atom below 2 angstr as 1, could be customized
             dist=sqrt(distx*distx+disty*disty+distz*distz)/3.7807
             write(std_out,*) iatom,iatom1,dist
             if (ab_mover%typat(iatom).ne.ab_mover%typat(iatom1)) then
               mcfac=-1
             else
               mcfac=1
             end if
             if (dist<1.0_dp)  dist=1.0_dp
             pot_perm(iatom)=pot_perm(iatom)+mcfac*(ab_mover%signperm)*1.0_dp&
&             /exp(log(dist)*6.0_dp)
           end if
         end do
       end do
       write(std_out,*) ' Perm_pot ',pot_perm(:)
!      write(message, '(a,10f12.5)' )' Perm_pot ',&
!      &         (pot_perm(iatom1),iatom1=1,ab_mover%natom)
!      call wrtout(ab_out,message,'COLL')
!      call wrtout(std_out,message,'COLL')

!      Find the two atoms, of different types, with the highest perm_pot
       max_perm(:)=-1.0d9
       do iatom=1,ab_mover%natom
         if (pot_perm(iatom) > max_perm(ab_mover%typat(iatom))) then
           max_perm(ab_mover%typat(iatom))=pot_perm(iatom)
           imax_perm(ab_mover%typat(iatom))=iatom
         end if
       end do

       if(zDEBUG)then
!        write(message, '(a,10f12.5)' )' max_Perm ',&
!        &      (max_perm(itypat),itypat=1,ab_mover%ntypat)
!        call wrtout(std_out,message,'COLL')
!        write(message, '(a,10i12)' )' imax_Perm ',&
!        &      (imax_perm(itypat),itypat=1,ab_mover%ntypat)
!        call wrtout(std_out,message,'COLL')
         write(std_out,*) 'NTYPAT',ab_mover%ntypat
         write(message, '(a,10f12.5)' )' max_Perm ',&
&         (max_perm(:))
         call wrtout(std_out,message,'COLL')
         write(message, '(a,10i12)' )' imax_Perm ',&
&         (imax_perm(:))
         call wrtout(std_out,message,'COLL')
       end if

!      Loop and keep the 2 largest values
       if (max_perm(1)>max_perm(2)) then
         maxp1=max_perm(1)
         maxp2=max_perm(2)
         iatom1=imax_perm(1)
         iatom2=imax_perm(2)
       else
         maxp1=max_perm(2)
         maxp2=max_perm(1)
         iatom1=imax_perm(2)
         iatom2=imax_perm(1)
       end if

       do itypat=3,ab_mover%ntypat
         if (max_perm(itypat)>maxp1) then
           maxp2=maxp1
           iatom2=iatom1
           maxp1=max_perm(itypat)
           iatom1=imax_perm(itypat)
         else if (max_perm(itypat)>maxp2) then
           maxp2=max_perm(itypat)
           iatom2=imax_perm(itypat)
         end if
       end do
       write(message, '(2(a,i5))' )' Will commute atom...',iatom1,'...of type ',&
&       ab_mover%typat(iatom1)
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
       write(message, '(2(a,i5))' )'         with atom...',iatom2,'...of type ',&
&       ab_mover%typat(iatom2)
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')

!      Commute the atoms positions
       distx=xcart(1,iatom1)
       disty=xcart(2,iatom1)
       distz=xcart(3,iatom1)
       xcart(1,iatom1)=xcart(1,iatom2)
       xcart(2,iatom1)=xcart(2,iatom2)
       xcart(3,iatom1)=xcart(3,iatom2)
       xcart(1,iatom2)=distx
       xcart(2,iatom2)=disty
       xcart(3,iatom2)=distz
!      Convert back to xred (reduced coordinates)
       call xcart2xred(ab_mover%natom,rprimd,xcart,xred)

     end if ! if(mod(itime,ab_mover%delayperm)==0)

   else
!    write(std_out,*) "I will jump to the end of cycle"
     jump_end_of_cycle=.TRUE.

   end if ! if (itime/=1 .and. ab_mover%delayperm/=0 .and. ab_mover%ntypat>2)
!  End of the commutation section

 end if ! if (icycle==1)

 if (allocated(imax_perm))   then
   ABI_DEALLOCATE(imax_perm)
 end if

!write(std_out,*) 'langevin 05',jump_end_of_cycle
!##########################################################
!### 05. Compute the next values (Only for extra cycles)

 if (icycle>1) then

!  write(std_out,*) "Entering internal cycle 2",icycle

!  If the energy computed is higher than the current
!  (etotal_temp) we have to discard the changes
!  and compute again

!  write(std_out,*) ch10
!  write(std_out,*) 'EVALUATION FORCES',etotal,hist%etot(abihist_findIndex(hist,-1))
!  write(std_out,*) ch10

!  This is the worst case (2 evaluations of SCFCV)
   ihist_prev = abihist_findIndex(hist,-1)
   if (etotal>hist%etot(ihist_prev).and.icycle==2) then

!    Discard the changes
     acell(:)   =hist%acell(:,ihist_prev)
     rprimd(:,:)=hist%rprimd(:,:,ihist_prev)
     xred(:,:)  =hist%xred(:,:,ihist_prev)
     fcart(:,:) =hist%fcart(:,:,ihist_prev)
     strten(:)  =hist%strten(:,ihist_prev)
     vel(:,:)   =hist%vel(:,:,ihist_prev)
     etotal     =hist%etot(ihist_prev)
     call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

!    distx=xcart(1,iatom1)
!    disty=xcart(2,iatom1)
!    distz=xcart(3,iatom1)
!    xcart(1,iatom1)=xcart(1,iatom2)
!    xcart(2,iatom1)=xcart(2,iatom2)
!    xcart(3,iatom1)=xcart(3,iatom2)
!    xcart(1,iatom2)=distx
!    xcart(2,iatom2)=disty
!    xcart(3,iatom2)=distz

!    Convert back to xred (reduced coordinates)
!    call xcart2xred(ab_mover%natom,rprimd,xcart,xred)
     write(message, '(a)' )' Commutation unsuccessful, recomputing the forces'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

!    This is the best case (only 1 evaluation of SCFCV)
   else

     write(message, '(a)')' Commutation successful ! Going on'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

!    Get rid of mean force on whole unit cell, but only if no generalized
!    constraints are in effect
!    call fcart2fred(fcart,fred_corrected,rprimd,ab_mover%natom)
!    if(ab_mover%nconeq==0)then
!      amass_tot=sum(ab_mover%amass(:))
!      do ii=1,3
!        if (ii/=3.or.ab_mover%jellslab==0) then
!          favg=sum(fred_corrected(ii,:))/dble(ab_mover%natom)
!          fred_corrected(ii,:)=fred_corrected(ii,:)-favg*ab_mover%amass(:)/amass_tot
!        end if
!      end do
!    end if

!    In thisc case we do not need to compute again SCFCV
!    We avoid the second iteration on ii
     jump_end_of_cycle=.TRUE.

   end if ! etotal > etotal_temp

 end if ! if (icycle>1)

!write(std_out,*) 'langevin 06',jump_end_of_cycle
!##########################################################
!### 06. Compute the next values (Only for the last cycle)

 if(jump_end_of_cycle)then
!  icycle=ncycle
!  write(std_out,*) 'This is the last cycle, avoid the others and continue'
   skipcycle=.TRUE.
 else
   skipcycle=.FALSE.
 end if

 if ((icycle==ncycle).OR.(skipcycle)) then

!  write(std_out,*) 'ENTERING THE FINAL PART',icycle,ncycle

!  Specific to Langevin dynamics
!  Initialize an array of random forces
!  No random force at itime=0
!  if (itime==0) then
   if (itime<0) then

     ran_force(:,:)=0.0_dp

   else

     do iatom=1,ab_mover%natom
!      sig_gauss is the std deviation of the random distribution
       sig_gauss=sqrt(2.0_dp*(ab_mover%friction)*ab_mover%amass(iatom)*ktemp)
!      write(std_out,*) 'sig_gauss=',sig_gauss
!      write(std_out,*) 'friction=',ab_mover%friction
!      write(std_out,*) 'ktemp=',ktemp
       do idim=1,3
         delxi=2.0_dp
         do while (delxi >= 1.0_dp)
           ran_num1=2.0_dp*uniformrandom(idum)-1.0_dp
           ran_num2=2.0_dp*uniformrandom(idum)-1.0_dp
           delxi=ran_num1*ran_num1+ran_num2*ran_num2
!          write(std_out,*) delxi,ran_num1,ran_num2
         end do
         ran_force(idim,iatom)=ran_num1*sqrt(-2.0_dp*log(delxi)/delxi)&
&         *sig_gauss/sqrt(ab_mover%dtion)
!        write(std_out,*) 'ran_force',ran_force(idim,iatom)
       end do
     end do

     if (zDEBUG) then
       write (std_out,*) '4. Different forces computed'
       write (std_out,*) 'RAN_FORCE'
       do kk=1,ab_mover%natom
         write (std_out,*) ran_force(:,kk)
       end do
     end if


     if(zDEBUG)then
!      The distribution should be gaussian
       delxi=0.0_dp
       do iatom=1,ab_mover%natom
         do idim=1,3
           delxi=delxi+(ran_force(idim,iatom)*ab_mover%dtion)**2
         end do
       end do
       delxi=delxi/(3.0_dp*ab_mover%natom)
       write(message, '(2(a,es22.14))' )' variance =',delxi,'  asked =',&
&       2.0_dp*(ab_mover%friction)*ab_mover%amass(2)*ktemp*ab_mover%dtion
       call wrtout(std_out,message,'COLL')
     end if
!    end if itime\=0

   end if

!  zDEBUG
!  write(message, '(a)' )' after initializing ran_force'
!  call wrtout(ab_out,message,'COLL')
!  call wrtout(std_out,message,'COLL')
!  ENDzDEBUG

   do iatom=1,ab_mover%natom
     do idim=1,3
       fcart_m(idim,iatom)=fcart(idim,iatom)/ab_mover%amass(iatom)
       ran_force(idim,iatom)=ran_force(idim,iatom)/ab_mover%amass(iatom)
     end do
   end do
   lang_force(:,:)=ran_force(:,:)-(ab_mover%friction)*vel(:,:)+fcart_m(:,:)

   if (zDEBUG) then
     write (std_out,*) '4. Different forces computed'
     write (std_out,*) 'FCART_M'
     do kk=1,ab_mover%natom
       write (std_out,*) fcart_m(:,kk)
     end do
     write (std_out,*) 'RAN_FORCE'
     do kk=1,ab_mover%natom
       write (std_out,*) ran_force(:,kk)
     end do
     write (std_out,*) 'LANG_FORCE'
     do kk=1,ab_mover%natom
       write (std_out,*) lang_force(:,kk)
     end do
   end if

!  zDEBUG
!  write(message, '(a)' )'before verlet'
!  call wrtout(ab_out,message,'COLL')
!  call wrtout(std_out,message,'COLL')
!  ENDzDEBUG

!  Compute next atomic coordinates using Verlet algorithm

!  Impose no change of acell, ucvol, rprim, and rprimd
   acell_next(:)=acell(:)
   ucvol_next=ucvol
   rprim_next(:,:)=rprim(:,:)
   rprimd_next(:,:)=rprimd(:,:)

!  Convert input xred (reduced coordinates) to xcart (cartesian)
   call xred2xcart(ab_mover%natom,rprimd,xcart,xred)
!  Uses the velocity
!
!  If an atom wants to cross the walls, velocity is reversed.
!
   do iatom=1,ab_mover%natom
     do idim=1,3
       delxi=xcart(idim,iatom)+ab_mover%dtion*vel(idim,iatom)+ &
&       0.5_dp*ab_mover%dtion*ab_mover%dtion*lang_force(idim,iatom)
       if ( (delxi > (rprimd(idim,idim)+(ab_mover%mdwall)) ) .or. &
&       (delxi < - (ab_mover%mdwall)                   )       ) then
         vel(idim,iatom)=-vel(idim,iatom)
         delxi=xcart(idim,iatom)+ab_mover%dtion*vel(idim,iatom)+ &
&         0.5_dp*ab_mover%dtion*ab_mover%dtion*lang_force(idim,iatom)
       end if
       xcart_next(idim,iatom)=delxi
     end do
   end do
   xcart(:,:)=xcart_next(:,:)
   if (zDEBUG) then
     write (std_out,*) '5. If an atom wants to cross the walls, velocity is reversed.'
     write (std_out,*) 'XCART'
     do kk=1,ab_mover%natom
       write (std_out,*) xcart(:,kk)
     end do
   end if

!  Convert back to xred_next (reduced coordinates)
   call xcart2xred(ab_mover%natom,rprimd,xcart_next,xred_next)
   call xcart2xred(ab_mover%natom,rprimd,xcart,xred)

   if (itime==1) then
!    no old forces are available at first step
!    Simple update of the velocity
!    first compute vel_nexthalf for next steps
     vel(:,:)=vel(:,:)+ab_mover%dtion*lang_force(:,:)
   else
!    case itime /= 0 normal verlet integration
     vel(:,:)=vel(:,:)+0.5_dp*ab_mover%dtion*(fcart_mold(:,:)+lang_force(:,:))
   end if
   if (zDEBUG) then
     write (std_out,*) '5. Change velocity with verlet'
     write (std_out,*) 'VEL'
     do kk=1,ab_mover%natom
       write (std_out,*) vel(:,kk)
     end do
   end if

!  Store 'current force' as 'old force'
   fcart_mold(:,:)=lang_force(:,:)

 end if ! if (icycle==ncycle)

!write(std_out,*) 'langevin 07',jump_end_of_cycle
!##########################################################
!### 07. Update the history with the prediction

!Increase indexes
 hist%ihist = abihist_findIndex(hist,+1)

 call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 hist%vel(:,:,hist%ihist)=vel(:,:)
 hist%time(hist%ihist)=real(itime,kind=dp)*ab_mover%dtion

 if (ab_mover%delayperm==0 .or. ab_mover%ntypat<=2) ncycle=1
 if(itime==ntime-1) ncycle=1

end subroutine pred_langevin
!!***

end module m_pred_langevin
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/pred_isothermal
!! NAME
!! pred_isothermal
!!
!! FUNCTION
!! Ionmov predictors (13) Isothermal integrator
!!
!! IONMOV 13:
!! Reversible integrator of Martyna at al.
!! The equation of motion of the ions in contact with a thermostat
!! and a barostat are solved with the algorithm proposed by Martyna,
!! Tuckermann Tobias and Klein [Mol. Phys., 1996, p. 1117].
!! Related parameters : the time step (dtion),
!! the initial temperature mdtemp(1), the final temperature mdtemp(2),
!! the number of thermostats (nnos), and the masses of thermostats (qmass).
!! If optcell=1 or 2, the mass of the barostat (bmass)
!! must be given in addition.
!!
!! There are three sub cases according to the value of optcell
!! optcell=0: isothermal
!! optcell=1: homogeneous cell fluctuations
!! optcell=2: full cell fluctuation in addition to temperature
!!            control.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, JCC, SE)
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
!! ntime  : Maximal number of iterations
!! zDEBUG : if true print some debugging information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces
!!                               acell, rprimd, stresses
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      dsyev,hist2var,isopress,isostress,isotemp,metric,mkrdim,var2hist,wrtout
!!      xcart2xred,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pred_isothermal(ab_mover,hist,itime,mttk_vars,ntime,zDEBUG,iexit)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_abimover
 use m_abihist
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pred_isothermal'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
 use interfaces_41_geometry
 use interfaces_45_geomoptim, except_this_one => pred_isothermal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(abimover),intent(in)       :: ab_mover
 type(abihist),intent(inout) :: hist
 type(mttk_type),intent(inout) :: mttk_vars
 integer,intent(in) :: itime
 integer,intent(in) :: ntime
 integer,intent(in) :: iexit
 logical,intent(in) :: zDEBUG

!Local variables-------------------------------
!scalars
 integer  :: ii,kk,iatom,idim,idum=5,ierr
 integer,parameter :: lwork=8
 real(dp) :: ucvol,ucvol0,ucvol_next,mttk_aloc,mttk_aloc2,mttk_bloc,ekin
 real(dp) :: massvol=0
 real(dp),parameter :: esh2=one/six,esh4=esh2/20._dp,esh6=esh4/42._dp
 real(dp),parameter :: esh8=esh6/72._dp,nosetol=tol10,v2tol=tol8
 real(dp) :: etotal,rescale_vel,polysh,s1,s2,sigma2,v2gauss,vtest
 real(dp),save :: ktemp,vlogv
 character(len=5000) :: message
!arrays
 real(dp),allocatable,save :: fcart_m(:,:),vel_nexthalf(:,:)

 real(dp) :: mttk_alc(3),mttk_alc2(3),mttk_blc(3),mttk_psh(3)
 real(dp) :: mttk_tv(3,3),mttk_vt(3,3),mttk_ubox(3,3)
 real(dp) :: mttk_uu(3),mttk_uv(3),mttk_veig(3)
 real(dp) :: acell(3),acell0(3),acell_next(3),favg_(3)
 real(dp) :: rprimd(3,3),rprimd0(3,3),rprim(3,3),rprimd_next(3,3),rprim_next(3,3)
 real(dp) :: gprimd(3,3)
 real(dp) :: gmet(3,3)
 real(dp) :: rmet(3,3)
 real(dp) :: fcart(3,ab_mover%natom)
!real(dp) :: fred_corrected(3,ab_mover%natom)
 real(dp) :: xcart(3,ab_mover%natom),xcart_next(3,ab_mover%natom)
 real(dp) :: xred(3,ab_mover%natom),xred_next(3,ab_mover%natom)
 real(dp) :: vel(3,ab_mover%natom)
 real(dp) :: strten(6),work(lwork)

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
   return
 end if

!write(std_out,*) 'isothermal 01'
!##########################################################
!### 01. Debugging and Verbose

 if(zDEBUG)then
   write(std_out,'(a,3a,41a,36a)') ch10,('-',kk=1,3),&
&   'Debugging and Verbose for pred_isothermal',('-',kk=1,36)
   write(std_out,*) 'ionmov: ',13
   write(std_out,*) 'itime:  ',itime
 end if

!write(std_out,*) 'isothermal 01'
!##########################################################
!### 01. Allocate the vectors vin, vout and hessian matrix
!###     These arrays could be allocated from a previus
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

!write(std_out,*) 'isothermal 02'
!##########################################################
!### 02. Obtain the present values from the history

 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)

 fcart(:,:)=hist%fcart(:,:,hist%ihist)
 strten(:) =hist%strten(:,hist%ihist)
 vel(:,:)  =hist%vel(:,:,hist%ihist)
 etotal    =hist%etot(hist%ihist)

 do ii=1,3
   rprim(ii,1:3)=rprimd(ii,1:3)/acell(1:3)
 end do
 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

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

!Save initial values
 acell0(:)=acell(:)
 rprimd0(:,:)=rprimd(:,:)
 ucvol0=ucvol

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

!write(std_out,*) 'isothermal 03'
!##########################################################
!### 05. Seconde half velocity step

 if (itime>1) then

!  Next Half velocity step
   do idim=1,3
     fcart_m(idim,:)=fcart(idim,:)/ab_mover%amass(:)
   end do
   vel(:,:)=vel_nexthalf(:,:)+ab_mover%dtion/two*fcart_m(:,:)

   if (ab_mover%optcell==0) then
!    Update Thermostat variables and velocity
     call isotemp(ab_mover%amass,ab_mover%dtion,ekin,ab_mover%iatfix,&
&     ktemp,mttk_vars,ab_mover%natom,ab_mover%nnos,ab_mover%qmass,vel)
   else if (ab_mover%optcell==1) then
!    Update Thermostat variables and velocity
     call isopress(ab_mover%amass,ab_mover%bmass,ab_mover%dtion,ekin,ab_mover%iatfix,&
&     ktemp,mttk_vars,ab_mover%natom,ab_mover%nnos,ab_mover%qmass,&
&     strten,ab_mover%strtarget,ucvol,vel,vlogv)
   else if (ab_mover%optcell==2) then
!    Next half step for extended variables
     call isostress(ab_mover%amass,ab_mover%bmass,ab_mover%dtion,ekin,ab_mover%iatfix,&
&     ktemp,mttk_vars,ab_mover%natom,ab_mover%nnos,&
&     ab_mover%qmass,strten,ab_mover%strtarget,ucvol,vel)
   end if

   if(itime==2) massvol=ekin+etotal

   if (ab_mover%optcell==2) then
!    Evolution of cell and volume
     acell_next(:)=acell(:)
     ucvol_next=ucvol
   end if

   call mkrdim(acell,rprim,rprimd)
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

   if(zDEBUG)then
     write(std_out,*) 'Second half velocity step'
     write(std_out,*) 'Cell parameters:'
     write(std_out,*) 'rprimd:'
     do kk=1,3
       write(std_out,*) rprimd(:,kk)
     end do
     write(std_out,*) 'rprim:'
     do kk=1,3
       write(std_out,*) rprim(:,kk)
     end do
     write(std_out,*) 'acell:'
     write(std_out,*) acell(:)
     write(std_out,*) 'Conserved energy:',(ekin+etotal)-massvol,ekin,etotal
     write(std_out,*) 'Volume of unitqry cell (ucvol):',ucvol
   end if

 end if ! if (itime>1)

!write(std_out,*) 'isothermal 04'
!##########################################################
!### 03. Compute the next values

!The temperature is linear between initial and final values
!It is here converted from Kelvin to Hartree (kb_HaK)
 ktemp=(ab_mover%mdtemp(1)+((ab_mover%mdtemp(2)-ab_mover%mdtemp(1))/dble(ntime-1))*(itime-1))*kb_HaK

 if(zDEBUG)then
   write(std_out,*) 'Temperature in Kelvin (ktemp):',ktemp
   write(std_out,*) 'Initial temp (mdtemp(1)):',ab_mover%mdtemp(1)
   write(std_out,*) 'Final temp (mdtemp(2)):',ab_mover%mdtemp(2)
   write(std_out,*) 'Delay for atom permutation (delayperm)',ab_mover%delayperm
   write(std_out,*) 'nnos:', ab_mover%nnos
   write(std_out,*) 'qmass', ab_mover%qmass(:)
   write(std_out,*) 'bmass',ab_mover%bmass
 end if

 if(itime==1) then
   mttk_vars%glogs(:)=zero; mttk_vars%vlogs(:)=zero; mttk_vars%xlogs(:)=zero
   mttk_vars%vboxg(:,:)=zero
   vlogv=zero
!  v2gauss is twice the kinetic energy
   v2gauss=0.0_dp
   do iatom=1,ab_mover%natom
     do idim=1,3
       v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*ab_mover%amass(iatom)
     end do
   end do

!  If there is no kinetic energy
   if (v2gauss<=v2tol.and.itime==1) then
!    Maxwell-Boltzman distribution
     v2gauss=zero
     vtest=zero
     do iatom=1,ab_mover%natom
       do idim=1,3
         vel(idim,iatom)=sqrt(kb_HaK*ab_mover%mdtemp(1)/ab_mover%amass(iatom))*cos(two_pi*uniformrandom(idum))
         vel(idim,iatom)=vel(idim,iatom)*sqrt(-2._dp*log(uniformrandom(idum)))
       end do
     end do

!    Get rid of center-of-mass velocity
     s1=sum(ab_mover%amass(:))
     do idim=1,3
       s2=sum(ab_mover%amass(:)*vel(idim,:))
       vel(idim,:)=vel(idim,:)-s2/s1
     end do

!    Recompute v2gauss
     do iatom=1,ab_mover%natom
       do idim=1,3
         v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*ab_mover%amass(iatom)
         vtest=vtest+vel(idim,iatom)/(3._dp*ab_mover%natom)
       end do
     end do

!    Now rescale the velocities to give the exact temperature
     rescale_vel=sqrt(3._dp*ab_mover%natom*kb_HaK*ab_mover%mdtemp(1)/v2gauss)
     vel(:,:)=vel(:,:)*rescale_vel

!    Recompute v2gauss with the rescaled velocities
     v2gauss=zero
     do iatom=1,ab_mover%natom
       do idim=1,3
         v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*ab_mover%amass(iatom)
       end do
     end do

!    Compute the variance and print
     sigma2=(v2gauss/(3._dp*ab_mover%natom)-ab_mover%amass(1)*vtest**2)/kb_HaK

     if (zDEBUG)then
       write(message, '(a)' )&
&       ' Rescaling or initializing velocities to initial temperature'
       call wrtout(std_out,message,'COLL')
       write(message, '(a,d12.5,a,D12.5)' )&
&       ' --- Scaling factor :',rescale_vel,' Asked T (K) ',ab_mover%mdtemp(1)
       call wrtout(std_out,message,'COLL')
       write(message, '(a,d12.5,a,D12.5)' )&
&       ' --- Effective temperature',v2gauss/(3*ab_mover%natom*kb_HaK),' From variance', sigma2
       call wrtout(std_out,message,'COLL')
     end if

   end if !(v2gauss<=v2tol.and.itime==1)
 end if !(itime==1)

!XG070613 : Do not take away the following line , seems needed for the pathscale compiler

 if (zDEBUG) write(std_out,*) 'vboxg',mttk_vars%vboxg(:,:)


!write(std_out,*) 'isothermal 05'
!##########################################################
!### 03. First half velocity step

!write(std_out,*) 'FIRST HALF VELOCITY STEP',ucvol
!write(std_out,*) 'OPTCELL option selected:',ab_mover%optcell
!write(std_out,*) 'RPRIMD'
!do kk=1,3
!write(std_out,*) rprimd(:,kk)
!end do
!write(std_out,*) 'RPRIM'
!do kk=1,3
!write(std_out,*) rprim(:,kk)
!end do
!write(std_out,*) 'ACELL'
!write(std_out,*) acell(:)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% BEGIN sub case optcell=0 Isothermal Ensemble
!%%%
 if(ab_mover%optcell==0) then
!  There is no evolution of cell
   acell_next(:)=acell(:)
   ucvol_next=ucvol
   rprim_next(:,:)=rprim(:,:)
   rprimd_next(:,:)=rprimd(:,:)
!  Update Thermostat variables and scale velocitie
   call isotemp(ab_mover%amass,ab_mover%dtion,ekin,ab_mover%iatfix,&
&   ktemp,mttk_vars,ab_mover%natom,ab_mover%nnos,ab_mover%qmass,vel)

!  Half velocity step
   do idim=1,3
     fcart_m(idim,:)=fcart(idim,:)/ab_mover%amass(:)
   end do
   vel_nexthalf(:,:)=vel(:,:)+ab_mover%dtion/two*fcart_m(:,:)
!  New positions
!  Convert input xred (reduced coordinates) to xcart (cartesian)
   call xred2xcart(ab_mover%natom,rprimd,xcart,xred)
   xcart_next(:,:)=xcart(:,:)+vel_nexthalf(:,:)*ab_mover%dtion
!  Convert back to xred (reduced coordinates)
   call xcart2xred(ab_mover%natom,rprimd,xcart_next,xred_next)
!  %%%
!  %%% END sub case optcell=0 Isothermal Ensemble
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%% BEGIN sub case optcell=1 Isothermal-Isenthalpic
!  %%%       Ensemble (homogeneous cell deformation)
!  %%%
 else if (ab_mover%optcell==1) then
!  Only homogeneous evolution of cell
!  Evolution of cell we keep rprim constant
   rprim_next(:,:)=rprim(:,:)
!  Update Thermostat variables and velocity
   call isopress(ab_mover%amass,ab_mover%bmass,ab_mover%dtion,ekin,ab_mover%iatfix,&
&   ktemp,mttk_vars,ab_mover%natom,ab_mover%nnos,ab_mover%qmass,&
&   strten,ab_mover%strtarget,ucvol,vel,vlogv)

!  Half velocity step
   do idim=1,3
     fcart_m(idim,:)=fcart(idim,:)/ab_mover%amass(:)
   end do
   vel_nexthalf(:,:)=vel(:,:)+ab_mover%dtion/two*fcart_m(:,:)
!  New positions
   mttk_aloc=exp(ab_mover%dtion/two*vlogv)
   mttk_aloc2=(vlogv*ab_mover%dtion/two)**2
   polysh=(((esh8*mttk_aloc2+esh6)*mttk_aloc2+esh4)*mttk_aloc2+esh2)*mttk_aloc2+one
   mttk_bloc=mttk_aloc*polysh*ab_mover%dtion
!  Convert input xred (reduced coordinates) to xcart (cartesian)
   call xred2xcart(ab_mover%natom,rprimd,xcart,xred)
   xcart_next(:,:)=xcart(:,:)*mttk_aloc**2+vel_nexthalf(:,:)*mttk_bloc
!  Update the volume and related quantities
   acell_next(:)=acell(:)*exp(ab_mover%dtion*vlogv)
!  ucvol=ucvol*exp(ab_mover%dtion*vlogv)
   call mkrdim(acell_next,rprim,rprimd_next)
   call metric(gmet,gprimd,-1,rmet,rprimd_next,ucvol_next)
!  Convert back to xred (reduced coordinates)
   call xcart2xred(ab_mover%natom,rprimd_next,xcart_next,xred_next)
!  Computation of the forces for the new positions
!  Compute LDA forces (big loop)

!  COMMENTED
!  This should be in mover.F90

!  !      If metric has changed since the initialization, update the Ylm's
!  if (ab_mover%optcell/=0.and.psps%useylm==1.and.itime>1)then
!  call status(0,dtfil%filstat,iexit,level,'call initylmg ')
!  option=0;if (ab_mover%iscf>0) option=1
!  call initylmg(gprimd,kg,ab_mover%kptns,ab_mover%mkmem,mpi_enreg,psps%mpsang,ab_mover%mpw,ab_mover%nband,ab_mover%nkpt,&
!  &         npwarr,ab_mover%nsppol,option,rprimd_next,ylm,ylmgr)
!  end if


!  %%%
!  %%% END sub case optcell=1 Isothermal-Isenthalpic
!  %%%     Ensemble (homogeneous cell deformation)
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%% BEGIN sub case optcell=2 Isothermal-Isenthalpic
!  %%%       Ensemble (full cell deformation)
!  %%%
 else if (ab_mover%optcell==2) then
   acell_next=acell
!  Fisrt half step for extended variables
   call isostress(ab_mover%amass,ab_mover%bmass,ab_mover%dtion,ekin,ab_mover%iatfix,&
&   ktemp,mttk_vars,ab_mover%natom,ab_mover%nnos,&
&   ab_mover%qmass,strten,ab_mover%strtarget,ucvol,vel)
!  Half velocity step
   do idim=1,3
     fcart_m(idim,:)=fcart(idim,:)/ab_mover%amass(:)
   end do
   vel_nexthalf(:,:)=vel(:,:)+ab_mover%dtion/two*fcart_m(:,:)
!  Convert input xred (reduced coordinates) to xcart (cartesian)
   call xred2xcart(ab_mover%natom,rprimd,xcart,xred)
!  New positions
   mttk_vt(:,:)=mttk_vars%vboxg(:,:)
   call dsyev('V','U',3,mttk_vt,3,mttk_veig,work,lwork,ierr)
   mttk_tv(:,:)=transpose(mttk_vt)
   mttk_alc(:)=exp(ab_mover%dtion/two*mttk_veig(:))
   mttk_alc2(:)=(mttk_veig(:)*ab_mover%dtion/two)**2
   mttk_psh(:)=(((esh8*mttk_alc2(:)+esh6)*mttk_alc2(:)+esh4)*mttk_alc2(:)+esh2)*mttk_alc2(:)+one
   mttk_blc(:)=mttk_alc(:)*mttk_psh(:)*ab_mover%dtion
!  Update the positions
   do iatom=1,ab_mover%natom
     mttk_uu(:)=matmul(mttk_tv,xcart(:,iatom))
     mttk_uv(:)=matmul(mttk_tv,vel_nexthalf(:,iatom))
     mttk_uu(:)=mttk_uu(:)*mttk_alc(:)**2+mttk_uv(:)*mttk_blc(:)
     xcart_next(:,iatom)=matmul(mttk_vt,mttk_uu)
   end do
!  Update the box (rprimd and rprim)
   mttk_ubox(:,:)=matmul(mttk_tv,rprimd)
   do idim=1,3
     mttk_ubox(:,idim)=mttk_ubox(:,idim)*mttk_alc(:)**2
   end do
   rprimd_next(:,:)=matmul(mttk_vt,mttk_ubox)
   do idim=1,3
     rprim_next(idim,:)=rprimd_next(idim,:)/acell(:)
   end do
!  Update the volume
   call metric(gmet,gprimd,-1,rmet,rprimd_next,ucvol)
!  Convert back to xred (reduced coordinates)
   call xcart2xred(ab_mover%natom,rprimd_next,xcart_next,xred_next)
!  Computation of the forces for the new positions

!  COMMENTED
!  This should be in mover.F90

!  !      If metric has changed since the initialization, update the Ylm's
!  if (ab_mover%optcell/=0.and.psps%useylm==1.and.itime>1)then
!  call status(0,dtfil%filstat,iexit,level,'call initylmg ')
!  option=0;if (ab_mover%iscf>0) option=1
!  call initylmg(gprimd,kg,ab_mover%kptns,ab_mover%mkmem,mpi_enreg,psps%mpsang,ab_mover%mpw,ab_mover%nband,ab_mover%nkpt,&
!  &         npwarr,ab_mover%nsppol,option,rprimd_next,ylm,ylmgr)
!  end if

!  %%%
!  %%% END sub case optcell=2 Isothermal-Isenthalpic
!  %%%     Ensemble (full cell deformation)
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 else
   write(message, '(a,i12,a,a)' )&
&   '  Disallowed value for optcell=',ab_mover%optcell,ch10,&
&   '  Allowed values with ionmov==13 : 0 to 2.'
   MSG_BUG(message)
 end if

!write(std_out,*) 'OLD PARAMETERS'
!write(std_out,*) 'RPRIMD'
!do kk=1,3
!write(std_out,*) rprimd(:,kk)
!end do
!write(std_out,*) 'RPRIM'
!do kk=1,3
!write(std_out,*) rprim(:,kk)
!end do
!write(std_out,*) 'ACELL'
!write(std_out,*) acell(:)

!write(std_out,*) 'NEXT PARAMETERS'
!write(std_out,*) 'RPRIMD'
!do kk=1,3
!write(std_out,*) rprimd_next(:,kk)
!end do
!write(std_out,*) 'RPRIM'
!do kk=1,3
!write(std_out,*) rprim_next(:,kk)
!end do
!write(std_out,*) 'ACELL'
!write(std_out,*) acell_next(:)


!Those are the values store into the history
 rprim=rprim_next
 rprimd=rprimd_next
 xred=xred_next
 xcart=xcart_next
 acell=acell_next

!write(std_out,*) 'isothermal 06'
!##########################################################
 !### 06. Update the history with the prediction

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

end subroutine pred_isothermal
!!***

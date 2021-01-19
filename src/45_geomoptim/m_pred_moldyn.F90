!!****m* ABINIT/m_pred_moldyn
!! NAME
!!  m_pred_moldyn
!!
!! FUNCTION
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

module m_pred_moldyn

 use defs_basis
 use m_abicore
 use m_abimover
 use m_abihist

 use m_geometry,  only : xcart2xred, xred2xcart
 use m_predtk,    only : fdtion

 implicit none

 private
!!***

 public :: pred_moldyn
!!***

contains
!!***

!!****f* ABINIT/pred_moldyn
!! NAME
!! pred_moldyn
!!
!! FUNCTION
!! Ionmov predictor (1) Molecular dynamics
!!
!! Molecular dynamics, with or without viscous damping
!! This function should be called after the call to scfcv
!! Updates positions, velocities and forces
!!
!! INPUTS
!! ab_mover<type abimover>=Subset of dtset only related with
!!          |                 movement of ions and acell, contains:
!!          | dtion:  Time step
!!          ! natom:  Number of atoms
!!          | vis:    viscosity
!!          | iatfix: Index of atoms and directions fixed
!!          | amass:  Mass of ions
!! icycle: Index of the internal cycle inside a time step (itime)
!! itime: Index of time iteration
!! zDEBUG : if true print some debugging information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist<type abihist>=Historical record of positions, forces,
!!                               stresses, cell and energies,
!!
!! ncycle: Number of cycles of a particular time step
!!
!! NOTES
!! * This routine is a predictor, it only produces new positions
!!   to be computed in the next iteration, this routine should
!!   produce not output at all
!! * ncycle changes from 4 for the first iteration (itime==1) to 1 for (itime>1)
!! * The arrays vec_tmp1 and vec_tmp2 are triky, they are use with
!!   different meanings, during the initialization they contains
!!   working positions and velocities that acumulated produce the
!!   first positions of itime=1, for itime>1 they will contain
!!   positions in 2 previous steps, those values are different
!!   from the values store in the history, thats the reason why
!!   we cannot simply use hist%xred to obtain those positions.
!!
!! PARENTS
!!      m_precpred_1geo
!!
!! CHILDREN
!!      hist2var,var2hist,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine pred_moldyn(ab_mover,hist,icycle,itime,ncycle,ntime,zDEBUG,iexit)

!Arguments ------------------------------------
!scalars
type(abimover),intent(in)       :: ab_mover
type(abihist),intent(inout),target :: hist
integer,intent(in)    :: icycle
integer,intent(inout) :: ncycle
integer,intent(in)    :: itime
integer,intent(in)    :: ntime
integer,intent(in)    :: iexit
logical,intent(in)    :: zDEBUG

!Local variables-------------------------------
!scalars
integer  :: kk,jj,ihist,ihist_next,ihist_prev,ihist_prev2
integer  :: ihist_prev4,ihist_prev5
real(dp) :: aa,alfa,bb,cc,x0,xm,em,vis,dx,dv
real(dp) :: fcart,fprev,fprev2
real(dp) :: xc
real(dp) :: vel,vnow,xnow,vprev
real(dp),save :: hh,time
!arrays
real(dp) :: acell(3),rprimd(3,3)
real(dp) :: xred(3,ab_mover%natom)
real(dp),allocatable :: xcart(:,:),xcart_prev(:,:)
real(dp),save,allocatable :: vec_tmp1(:,:)
real(dp),save,allocatable :: vec_tmp2(:,:)
real(dp), ABI_CONTIGUOUS pointer :: vel_cur(:,:),vel_next(:,:)
real(dp),pointer :: fcart_cur(:,:),fcart_prev(:,:),fcart_prev2(:,:)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 if(iexit/=0)then
   if(allocated(vec_tmp1))  then
     ABI_DEALLOCATE(vec_tmp1)
   end if
   if(allocated(vec_tmp2))  then
     ABI_DEALLOCATE(vec_tmp2)
   end if
   return
 end if

 vis= ab_mover%vis
!Just to avoid warnings of uninitialized variables
 fprev=0.0_dp
 fprev=0.0_dp
 fprev2=0.0_dp
 vnow=0.0_dp
 vprev=0.0_dp
 xnow=0.0_dp

!Those arrays contains intermediary results used with
!different meanings during the different time steps
!We need to preserv the allocation status, this is the
!reason to be 'SAVE'
 if (itime==1.and.icycle==1)then
   if(allocated(vec_tmp1))  then
     ABI_DEALLOCATE(vec_tmp1)
   end if
   if(allocated(vec_tmp2))  then
     ABI_DEALLOCATE(vec_tmp2)
   end if
   ABI_ALLOCATE(vec_tmp1,(3,ab_mover%natom))
   ABI_ALLOCATE(vec_tmp2,(3,ab_mover%natom))
 end if

!write(std_out,*) '00'
!##########################################################
!### 00. Copy from the history to the variables

 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)

 ABI_ALLOCATE(xcart,(3,ab_mover%natom))
 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

 if (itime==1.or.itime==2)then
   ABI_ALLOCATE(xcart_prev,(3,ab_mover%natom))
   call xred2xcart(ab_mover%natom,rprimd,xcart_prev,hist%xred(:,:,1))
 end if

 ihist = abihist_findIndex(hist, 0)
 ihist_prev  = abihist_findIndex(hist,-1)
 ihist_prev2 = abihist_findIndex(hist,-2)
 ihist_prev4 = abihist_findIndex(hist,-4)
 ihist_prev5 = abihist_findIndex(hist,-5)
 ihist_next  = abihist_findIndex(hist,+1)

 fcart_cur => hist%fcart(:,:,ihist)
 if (itime==2) fcart_prev  => hist%fcart(:,:,ihist_prev4)
 if (itime==3) fcart_prev2 => hist%fcart(:,:,ihist_prev5)
 if (itime >2.or. icycle>=2)fcart_prev  => hist%fcart(:,:,ihist_prev)
 if (itime >3.or. icycle>=3)fcart_prev2 => hist%fcart(:,:,ihist_prev2)

 vel_cur  => hist%vel(:,:,ihist)
 vel_next => hist%vel(:,:,ihist_next)

!write(std_out,*) '01'
!##########################################################
!### 01. Get or compute the time step dtion

 if (ab_mover%dtion>0)then
   hh = ab_mover%dtion
 else
   hh=fdtion(ab_mover,itime,xcart,fcart_cur,vel_cur)
 end if

!write(std_out,*) '02'
!##########################################################
!### 02. For all atoms and directions
 do kk=1,ab_mover%natom
   em=ab_mover%amass(kk)
   do jj=1,3

!    write(std_out,*) '03'
!    ##########################################################
!    ### 03. Filling other values from history (forces and vel)
     fcart=fcart_cur(jj,kk)
     xc=xcart(jj,kk)
     vel=hist%vel(jj,kk,1)

!    Previous values only after first iteration
     if (itime>=2.or.icycle>=2) then
       fprev=fcart_prev(jj,kk)
       vprev=hist%vel(jj,kk,hist%ihist)
     end if
     if (itime>=3.or.icycle>=3) then
       fprev2=fcart_prev2(jj,kk)
     end if

     if (itime==2)then
       vec_tmp1(jj,kk)=xcart_prev(jj,kk)
       vec_tmp2(jj,kk)=xcart(jj,kk)
     end if

!    write(std_out,*) '04'
!    ##########################################################
!    ### 04. Take first the atoms that are not allowed to move along
!    ###     this direction
!    ###     Warning : implemented in cartesian coordinates
     if (ab_mover%iatfix(jj,kk)==1) then
!      Their positions will be the same as xcart
       xnow=xcart(jj,kk)
!      Their velocities are zero
       vnow=0.0_dp
     else

!      write(std_out,*) '05'
!      ##########################################################
!      ### 05. Initialization (itime==1):
!      ###     4 calls to obtain the forces are neeeded
!      ###     The variables vec_tmp2 and vec_tmp1 from previous
!      ###     calls are used in the following ones.
       if(itime==1)then
         x0=xcart_prev(jj,kk)

!        Prepare the second cycle
         if(icycle==1)then
           dx=hh*vel
           dv=hh/em*(fcart-vis*vel)
           xnow=x0+.5_dp*dx
           vnow=vel+.5_dp*dv
           vec_tmp2(jj,kk)=xc+sixth*dx
           vec_tmp1(jj,kk)=vel+sixth*dv
         else if(icycle==2)then
           dx=hh*vprev
           dv=hh/em*(fcart-vis*vprev)
           xnow=x0+.5_dp*dx
           vnow=vel+.5_dp*dv
           vec_tmp2(jj,kk)=vec_tmp2(jj,kk)+third*dx
           vec_tmp1(jj,kk)=vec_tmp1(jj,kk)+third*dv
         else if(icycle==3)then
           dx=hh*vprev
           dv=hh/em*(fcart-vis*vprev)
           xnow=x0+dx
           vnow=vel+dv
           vec_tmp2(jj,kk)=vec_tmp2(jj,kk)+third*dx
           vec_tmp1(jj,kk)=vec_tmp1(jj,kk)+third*dv
         else if(icycle==4)then
           dx=hh*vprev
           dv=hh/em*(fcart-vis*vprev)
           xnow=vec_tmp2(jj,kk)+sixth*dx
           vnow=vec_tmp1(jj,kk)+sixth*dv
         end if
       else !(itime/=1)

!        write(std_out,*) '06'
!        ##########################################################
!        ### 06. Change positions and velocities
!        ###     These changes only applies for itime>2
         if (itime>2)then
!          Uses a corrector to have better value of xnow, and
!          derive vnow. Only update atoms position and
!          velocity along its allowed directions
           aa=fprev
           bb=(fcart-fprev2)/(2._dp*hh)
           cc=(fcart+fprev2-2._dp*fprev)/(2._dp*hh*hh)
           x0=vec_tmp2(jj,kk)
           xm=vec_tmp1(jj,kk)
           if(abs(vis)<=1.d-8)then
!            NON-DAMPED DYNAMICS (Post-Code)
             xnow=2._dp*x0-xm+hh**2/em/12._dp*&
&             (fprev2+10._dp*fprev+fcart)
             vnow=(bb*hh**2)/(3._dp*em)&
&             +1.5_dp*aa*hh/em+&
&             (5._dp/12._dp)*cc*hh**3/em&
&             +x0/hh-xm/hh
           else
!            DAMPED DYNAMICS (Post-Code)
             alfa=exp(-vis*hh/em)
             xnow=((-aa*hh*vis**2+0.5_dp*bb*hh**2*vis**2&
&             -third*cc*hh**3*vis**2+em*bb*hh*vis&
&             -em*cc*hh**2*vis-2._dp*em**2*cc*hh+x0*vis**3-xm*vis**3)*alfa&
&             +aa*hh*vis**2-em*bb*hh*vis+third*cc*hh**3*vis**2&
&             +2._dp*em**2*cc*hh+0.5D0*bb*hh**2*vis**2-em*cc*hh**2*vis+x0*vis**3)&
&             /vis**3
             vnow=(em*aa*vis**2*alfa-em*aa*vis**2+bb*hh*vis**2*em*alfa&
&             -bb*hh*vis**2*em+cc*hh**2*vis**2*em*alfa-cc*hh**2*vis**2*em&
&             -em**2*bb*vis*alfa+em**2*bb*vis-2._dp*em**2*cc*hh*vis*alfa+&
&             2._dp*em**2*cc*hh*vis+2._dp*em**3*cc*alfa-2._dp*em**3*cc+&
&             vis**3*alfa**2*aa*hh-0.5_dp*vis**3*alfa**2*bb*hh**2+&
&             third*vis**3*alfa**2*cc*hh**3-vis**2*&
&             alfa**2*em*bb*hh+vis**2*alfa**2*em*cc*hh**2+&
&             2._dp*vis*alfa**2*em**2*cc*hh-vis**4*alfa**2*x0+&
&             vis**4*alfa**2*xm)/vis**3/(alfa-1._dp)/em

           end if !if(abs(vis)<=1.d-8)

           xc=xnow
           vec_tmp1(jj,kk)=vec_tmp2(jj,kk)
           vec_tmp2(jj,kk)=xnow
         else
           vnow=vprev
         end if !if(itime>2)

!        write(std_out,*) '07'
!        ##########################################################
!        ### 07. Correct positions
!        ###     These changes only applies for itime>1

         if(abs(vis)<=1.d-8)then
!          NON-DAMPED DYNAMICS (Pre-Code)
!          If the viscosity is too small, the equations become
!          ill conditioned due to rounding error so do regular
!          Verlet predictor Numerov corrector.
           x0=vec_tmp2(jj,kk)
           xm=vec_tmp1(jj,kk)
           xnow=2._dp*x0-xm&
&           + hh**2/em*fcart
         else
!          DAMPED DYNAMICS (Pre-Code)
!          These equations come from solving
!          m*d2x/dt2+vis*dx/dt=a+b*t+c*t**2
!          analytically under the boundary conditions that
!          x(0)=x0 and x(-h)=xm, and the following is the
!          expression for x(h). a, b and c are determined
!          from our knowledge of the driving forces.
           aa=fcart
           bb=(fcart-fprev)/hh
           x0=vec_tmp2(jj,kk)
           xm=vec_tmp1(jj,kk)
           alfa=exp(-vis*hh/em)
           xnow=( (-aa*hh*vis**2 +0.5_dp*bb*hh**2*vis**2&
&           +em*bb*hh*vis +x0*vis**3 -xm*vis**3)*alfa&
&           +aa*hh*vis**2 -em*bb*hh*vis&
&           +0.5_dp*bb*hh**2*vis**2 +x0*vis**3)/vis**3
!          End of choice between initialisation, damped
!          dynamics and non-damped dynamics
         end if

       end if !if(itime==1)

     end if !if(ab_mover%iatfix(jj,kk)==1)

!    write(std_out,*) '08'
!    ##########################################################
!    ### 08. Update history

     xcart(jj,kk)=xnow
     vel_next(jj,kk)=vnow

!    write(std_out,*) '09'
!    ##########################################################
!    ### 09. End loops of atoms and directions
   end do ! jj=1,3
 end do ! kk=1,ab_mover%natom

!write(std_out,*) '10'
!##########################################################
!### 10. Filling history with the new values

 hist%ihist = abihist_findIndex(hist,+1)

 call xcart2xred(ab_mover%natom,rprimd,xcart,xred)
 call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)

!Change ncycle for itime>1
 if (icycle==4)  ncycle=1

 if (itime==1)then
   time=0.0_dp
   if (ab_mover%dtion<0)then
     write(std_out,*) 'Time=',time
   end if
 end if
 time=time+hh
 hist%time(hist%ihist)=time

 if (allocated(xcart)) then
   ABI_DEALLOCATE(xcart)
 end if
 if (allocated(xcart_prev)) then
   ABI_DEALLOCATE(xcart_prev)
 end if

 if (.false.) write(std_out,*) ntime

end subroutine pred_moldyn
!!***

end module m_pred_moldyn
!!***

!!****m* ABINIT/m_pred_steepdesc
!! NAME
!!  m_pred_steepdesc
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

module m_pred_steepdesc

 use defs_basis
 use m_abicore
 use m_abimover
 use m_abihist

 use m_geometry,       only : mkradim, mkrdim, xcart2xred, xred2xcart
 use m_predtk,         only : fdtion

 implicit none

 private
!!***

 public :: pred_steepdesc
!!***

contains
!!***

!!****f* ABINIT/pred_steepdesc
!! NAME
!! pred_steepdesc
!!
!! FUNCTION
!! Ionmov predictor (21) Steepest Descent Algorithm
!! The update of positions is given by the following equation:
!!
!! $$\Delta r_{n,i}=\lambda F_{n,i}$$
!!
!! r is the position of the 'n' ion along the 'i' direction
!! F is the force of the 'n' ion along the 'i' direction.
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
!! hist<type abihist>=Historical record of positions, forces, stresses, cell and energies.
!!
!! ncycle: Number of cycles of a particular time step
!!
!! NOTES
!! * This routine is a predictor, it only produces new positions
!!   to be computed in the next iteration, this routine should
!!   produce not output at all
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      hist2var,mkradim,mkrdim,var2hist,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine pred_steepdesc(ab_mover,forstr,hist,itime,zDEBUG,iexit)

implicit none

!Arguments ------------------------------------
!scalars
type(abimover),intent(in)       :: ab_mover
type(abihist),intent(inout),target :: hist
type(abiforstr),intent(in) :: forstr
integer,intent(in)    :: itime,iexit
logical,intent(in)    :: zDEBUG

!Local variables-------------------------------
!scalars
integer  :: kk,jj,ihist_prev
real(dp) :: em
real(dp) :: f_cart
real(dp) :: xc,str
real(dp) :: xnow,lambda
real(dp),save :: hh
!arrays
real(dp) :: acell(3),strten(6)
real(dp) :: rprim(3,3),rprimd(3,3)
real(dp) :: xred(3,ab_mover%natom),xcart(3,ab_mover%natom)
real(dp) :: residual(3,ab_mover%natom)
real(dp), ABI_CONTIGUOUS pointer :: fcart(:,:),vel(:,:)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 if(iexit/=0)then
   return
 end if

!write(std_out,*) '01'
!##########################################################
!### 01. Copy from the history to the variables
 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)

 do jj=1,3
   rprim(jj,1:3)=rprimd(jj,1:3)/acell(1:3)
 end do

 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)
 strten(:)=hist%strten(:,hist%ihist)
 fcart => hist%fcart(:,:,hist%ihist)
 vel => hist%vel(:,:,hist%ihist)

!Fill the residual with forces (No preconditioning)
!Or the preconditioned forces
 if (ab_mover%goprecon==0)then
   residual(:,:)=fcart(:,:)
 else
   residual(:,:)= forstr%fcart(:,:)
 end if

!write(std_out,*) '01'
!##########################################################
!### 02. Get or compute de time step dtion

 if (ab_mover%dtion>0)then
   hh = ab_mover%dtion
 else
   hh=fdtion(ab_mover,itime,xcart,fcart,vel)
 end if

 lambda=hh
 write(std_out,*) 'Lambda',lambda

!write(std_out,*) '02'
!##########################################################
!### 02. For all atoms and directions
 do kk=1,ab_mover%natom
!  Normally this is the mass of the atom
!  em=ab_mover%amass(kk)
!  But the steepest algorithm is using unitary mass
   em=1
   do jj=1,3

!    write(std_out,*) '03'
!    ##########################################################
!    ### 03. Filling other values from history (forces and vel)
     f_cart=residual(jj,kk)
     xc=xcart(jj,kk)
!    This lambda is taken from the kinematical equation
!    lambda=hh*hh/(2*em)

!    write(std_out,*) '04'
!    ##########################################################
!    ### 04. Take first the atoms that are not allowed to move along
!    ###     this direction
!    ###     Warning : implemented in cartesian coordinates
     if (ab_mover%iatfix(jj,kk)==1) then
!      Their positions will be the same as xcart
       xnow=xc
     else

!      This is the main expresion (1)
       xnow=xc+lambda*f_cart

     end if !if(ab_mover%iatfix(jj,kk)==1)

!    write(std_out,*) '05'
!    ##########################################################
!    ### 08. Update history

     xcart(jj,kk)=xnow

!    write(std_out,*) '06'
!    ##########################################################
!    ### 09. End loops of atoms and directions
   end do ! jj=1,3
 end do ! kk=1,ab_mover%natom

 if (ab_mover%optcell/=0)then

   if (ab_mover%optcell==1)then
     do jj=1,3
       acell(jj)=acell(jj)+lambda*strten(jj)
     end do ! jj=1,3
     call mkrdim(acell,rprim,rprimd)
   elseif (ab_mover%optcell==2)then
     do kk=1,3
       do jj=1,3
         if (jj==1 .and. kk==1) str=strten(1)
         if (jj==2 .and. kk==2) str=strten(2)
         if (jj==3 .and. kk==3) str=strten(3)
         if (jj==1 .and. kk==2) str=strten(6)
         if (jj==1 .and. kk==3) str=strten(5)
         if (jj==2 .and. kk==1) str=strten(6)
         if (jj==3 .and. kk==1) str=strten(5)
         if (jj==2 .and. kk==3) str=strten(4)
         if (jj==3 .and. kk==2) str=strten(4)
         rprimd(jj,kk)=rprimd(jj,kk)+lambda*str
       end do ! jj=1,3
     end do ! kk=1,3
     call mkradim(acell,rprim,rprimd)
   end if

 end if


!write(std_out,*) '08'
!##########################################################
!### 10. Filling history with the new values

!Increase indices
 hist%ihist = abihist_findIndex(hist,+1)

!Compute xred from xcart, and rprimd
 call xcart2xred(ab_mover%natom,rprimd,xcart,xred)

 call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 ihist_prev = abihist_findIndex(hist,-1)
 hist%vel(:,:,hist%ihist)=hist%vel(:,:,ihist_prev)

end subroutine pred_steepdesc
!!***

end module m_pred_steepdesc
!!***

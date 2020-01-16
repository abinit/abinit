!!****m* ABINIT/m_pred_srkhna14
!! NAME
!!  m_pred_srkna14
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

module m_pred_srkhna14

 use defs_basis
 use m_abicore
 use m_abimover
 use m_abihist

 use m_geometry,    only : xcart2xred, xred2xcart, metric

 implicit none

 private
!!***

 public :: pred_srkna14
!!***

contains
!!***

!!****f* ABINIT/pred_srkna14
!! NAME
!! pred_srkna14
!!
!! FUNCTION
!! Ionmov predictors (14) Srkna14 molecular dynamics
!!
!! IONMOV 14:
!! Simple molecular dynamics with a symplectic algorithm proposed
!! by S.Blanes and P.C.Moans, called SRKNa14 in Practical symplectic partitioned
!! Runge--Kutta and Runge--Kutta--Nystrom methods, Journal of Computational
!! and Applied Mathematics archive, volume 142,  issue 2  (May 2002), pages 313 - 330 [[cite:Blanes2002]].
!! of the kind first published by H. Yoshida, Construction of higher order symplectic
!! integrators, Physics Letters A, volume 150, number 5 to 7, pages 262 - 268 [[cite:Yoshida1990]]
!! This algorithm requires at least 14 evaluation of the forces (actually 15 are done
!! within Abinit) per time step. At this cost it usually gives much better
!! energy conservation than the verlet algorithm (ionmov 6) for a 30 times bigger
!! value of <a href="varrlx.html#dtion">dtion</a>. Notice that the potential
!! energy of the initial atomic configuration is never evaluated using this
!! algorithm.
!!
!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information needed by the preditor
!! itime  : Index of the present iteration
!! ntime  : Maximal number of iterations
!! icycle : Index of the present cycle
!! zDEBUG : if true print some debugging information
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces acell, rprimd, stresses
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      hist2var,metric,var2hist,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine pred_srkna14(ab_mover,hist,icycle,zDEBUG,iexit,skipcycle)

 implicit none

!Arguments ------------------------------------
!scalars
 type(abimover),intent(in)       :: ab_mover
 type(abihist),intent(inout) :: hist
 integer,intent(in)    :: icycle
 integer,intent(in)    :: iexit
 logical,intent(in)    :: zDEBUG
 logical,intent(out)   :: skipcycle

!Local variables-------------------------------
!scalars
 integer  :: ihist_prev,ii,jj,kk
 real(dp) :: ucvol,ucvol_next
 real(dp),parameter :: v2tol=tol8
 real(dp) :: etotal
 logical  :: jump_end_of_cycle=.FALSE.
! character(len=5000) :: message
!arrays
 real(dp),save :: aa(15),bb(15)
 real(dp) :: acell(3),acell_next(3)
 real(dp) :: rprimd(3,3),rprimd_next(3,3)
 real(dp) :: gprimd(3,3),gmet(3,3),rmet(3,3)
 real(dp) :: fcart(3,ab_mover%natom),fcart_m(3,ab_mover%natom)
!real(dp) :: fred_corrected(3,ab_mover%natom)
 real(dp) :: xcart(3,ab_mover%natom)
 real(dp) :: xred(3,ab_mover%natom)
 real(dp) :: vel(3,ab_mover%natom)
 real(dp) :: strten(6)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 if(iexit/=0)then
   return
 end if

 jump_end_of_cycle=.FALSE.
 fcart_m(:,:)=zero

!write(std_out,*) 'srkna14 03',jump_end_of_cycle
!##########################################################
!### 03. Obtain the present values from the history

 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

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

!write(std_out,*) 'srkna14 04',jump_end_of_cycle
!##########################################################
!### 04. Compute the next values (Only for the first cycle)

 if (icycle==1) then

   if(zDEBUG) then
     write(std_out,*) 'Entering only for first cycle'
   end if

   aa(1) =  0.0378593198406116_dp;
   aa(2) =  0.102635633102435_dp;
   aa(3) = -0.0258678882665587_dp;
   aa(4) =  0.314241403071447_dp;
   aa(5) = -0.130144459517415_dp;
   aa(6) =  0.106417700369543_dp;
   aa(7) = -0.00879424312851058_dp;
   aa(8) =  1._dp -&
&   2._dp*(aa(1)+aa(2)+aa(3)+aa(4)+aa(5)+aa(6)+aa(7));
   aa(9) =  aa(7);
   aa(10)=  aa(6);
   aa(11)=  aa(5);
   aa(12)=  aa(4);
   aa(13)=  aa(3);
   aa(14)=  aa(2);
   aa(15)=  aa(1);

   bb(1) =  0.0_dp
   bb(2) =  0.09171915262446165_dp;
   bb(3) =  0.183983170005006_dp;
   bb(4) = -0.05653436583288827_dp;
   bb(5) =  0.004914688774712854_dp;
   bb(6) =  0.143761127168358_dp;
   bb(7) =  0.328567693746804_dp;
   bb(8) =  0.5_dp - (bb(1)+bb(2)+bb(3)+bb(4)+bb(5)+bb(6)+bb(7));
   bb(9) =  0.5_dp - (bb(1)+bb(2)+bb(3)+bb(4)+bb(5)+bb(6)+bb(7));
   bb(10)=  bb(7);
   bb(11)=  bb(6);
   bb(12)=  bb(5);
   bb(13)=  bb(4);
   bb(14)=  bb(3);
   bb(15)=  bb(2);

   acell_next(:)=acell(:)
   ucvol_next=ucvol
   rprimd_next(:,:)=rprimd(:,:)

!  step 1 of 15

!  Convert input xred (reduced coordinates) to xcart (cartesian)
   call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

   vel(:,:) = vel(:,:) + bb(1) * ab_mover%dtion * fcart_m(:,:)

   do ii=1,3
     do jj=1,ab_mover%natom
       write(std_out,*) xcart(ii,jj), ab_mover%dtion, aa(1), vel(ii,jj)
       xcart(ii,jj) = xcart(ii,jj) + ab_mover%dtion * aa(1) * vel(ii,jj)
       write(std_out,*) xcart(ii,jj)
     end do
   end do

!  xcart(:,:) = xcart(:,:) + ab_mover%dtion * aa(1) * vel(:,:);

!  Convert back to xred (reduced coordinates)
   call xcart2xred(ab_mover%natom,rprimd,xcart,xred)

 end if ! if (icycle==1)

!write(std_out,*) 'srkna14 05',jump_end_of_cycle
!##########################################################
!### 05. Compute the next values (Only for extra cycles)

 if (icycle>1) then

   do ii=1,ab_mover%natom
     do jj=1,3
       fcart_m(jj,ii) = fcart(jj,ii)/ab_mover%amass(ii)
     end do
   end do

   if (icycle<16)then

!    Update of velocities and positions
     vel(:,:) = vel(:,:) + bb(icycle) * ab_mover%dtion * fcart_m(:,:)
     xcart(:,:) = xcart(:,:) +&
&     aa(icycle) * ab_mover%dtion * vel(:,:)
!    Convert xcart_next to xred_next (reduced coordinates)
!    for scfcv
     call xcart2xred(ab_mover%natom, rprimd, xcart,&
&     xred)

   end if ! (ii<16)

 end if ! if (icycle>1)

!write(std_out,*) 'srkna14 06',jump_end_of_cycle
!##########################################################
!### 06. Compute the next values (Only for the last cycle)

 if(jump_end_of_cycle)then
   skipcycle=.TRUE.
 else
   skipcycle=.FALSE.
 end if

!write(std_out,*) 'srkna14 07',jump_end_of_cycle
!##########################################################
!### 07. Update the history with the prediction

!Increase indexes
 hist%ihist = abihist_findIndex(hist,+1)

!Fill the history with the variables
!xred, acell, rprimd, vel
 call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 hist%vel(:,:,hist%ihist)=vel(:,:)
 ihist_prev = abihist_findIndex(hist,-1)
 hist%time(hist%ihist)=hist%time(ihist_prev)+ab_mover%dtion

end subroutine pred_srkna14
!!***

end module m_pred_srkhna14
!!***

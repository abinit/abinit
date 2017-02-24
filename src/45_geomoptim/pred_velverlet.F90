!{\src2tex{textfont=tt}}
!!****f* ABINIT/pred_velverlet
!! NAME
!!  pred_velverlet
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2017 ABINIT group (FIXME: add author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
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


subroutine pred_velverlet(ab_mover,hist,itime,ntime,zDEBUG,iexit)
    
 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_abimover
 use m_abihist

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pred_velverlet'
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(abimover),intent(in)   :: ab_mover
 type(abihist),intent(inout) :: hist
 integer,intent(in) :: itime
 integer,intent(in) :: ntime
 integer,intent(in) :: iexit
 logical,intent(in) :: zDEBUG

!Local variables-------------------------------

 integer  :: ii,jj
 real(dp) :: amass_tot,favg,etotal,epot,ekin,ekin_tmp
 real(dp) :: xcart(3,ab_mover%natom),xcart_next(3,ab_mover%natom)
 real(dp) :: xred(3,ab_mover%natom),xred_next(3,ab_mover%natom)
 real(dp) :: vel(3,ab_mover%natom),vel_nexthalf(3,ab_mover%natom)
 real(dp) :: fcart(3,ab_mover%natom),fred(3,ab_mover%natom)
 real(dp) :: fred_corrected(3,ab_mover%natom),fcart_corrected(3,ab_mover%natom)

 real(dp) :: acell(3)
 real(dp) :: rprimd(3,3),rprim(3,3)
!real(dp),allocatable,save :: vin_prev(:)
! *************************************************************************

 DBG_ENTER("COLL")
 
! if (option/=1 .and. option/=2 ) then
!   write(msg,'(3a,i0)')&
!&   'The argument option should be 1 or 2,',ch10,&
!&   'however, option=',option
!   MSG_BUG(msg)
! end if
!
! if (sizein<1) then
!   write(msg,'(3a,i0)')&
!&   'The argument sizein should be a positive number,',ch10,&
!&   'however, sizein=',sizein
!   MSG_ERROR(msg)
! end if

 DBG_EXIT("COLL")

if(iexit/=0)then
  return
end if

 ! Start preparation for velocity verlet

 call hist2var(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,zDEBUG)

 fcart(:,:) = hist%histXF(:,:,3,hist%ihist)
 fred(:,:)  = hist%histXF(:,:,4,hist%ihist)
 vel(:,:)   = hist%histV(:,:,hist%ihist)               
 epot       = hist%histE(hist%ihist)                   ! electronic sub-system energy
 ekin       = hist%histE(hist%ihist)                   ! kinetic energy


 if(zDEBUG)then
   write (std_out,*) 'velverlet step',itime
   write (std_out,*) 'fcart:'
   do ii=1,3
     write (std_out,*) fcart(ii,:)
   end do
   write (std_out,*) 'fred:'
   do ii=1,3
     write (std_out,*) fred(ii,:)
   end do
   write (std_out,*) 'xcart:'
   do ii=1,3
     write (std_out,*) xcart(ii,:)
   end do
   write (std_out,*) 'vel:'
   do ii=1,3
     write (std_out,*) vel(ii,:)
   end do
 end if


 do ii=1,ab_mover%natom ! propagate velocities half time step forward
   do jj=1,3
     vel(jj,ii) = vel(jj,ii) + 0.5_dp * ab_mover%dtion*fcart(jj,ii)/ab_mover%amass(ii)
   enddo
 enddo

 if(itime/=1)then
    ! Here, the time moments of velocities and coordinates coincide
    ! Compute the kinetic energy
    ekin_tmp=0.0
    do ii=1,ab_mover%natom
      do jj=1,3
        ekin_tmp=ekin_tmp+0.5_dp*ab_mover%amass(ii)*vel(jj,ii)**2
      end do
    end do
    !and print it out for debug (full energy, i.e. kinetic + potential should not drift)
    if(zDEBUG)then
          write (std_out,*) 'epot,ekin:', epot,ekin_tmp
    endif 

    if(itime/=ntime)then
     do ii=1,ab_mover%natom ! propagate velocities another half time step forward 
        do jj=1,3           ! "if" condition established a shift of half time step between velocity and coordinates time grids
           vel(jj,ii) = vel(jj,ii) + 0.5_dp * ab_mover%dtion*fcart(jj,ii)/ab_mover%amass(ii)
        enddo
     enddo
    endif
 endif


 if(itime/=ntime)then
  do ii=1,ab_mover%natom ! propagate coordinates one time step forward
    do jj=1,3
      xcart(jj,ii) = xcart(jj,ii) + ab_mover%dtion*vel(jj,ii)
    enddo
  enddo
 endif

! convert new xcart to xred to set correct output values 
! update the history with the new coordinates, velocities, etc.

!Increase indexes
 hist%ihist=hist%ihist+1

 call xcart2xred(ab_mover%natom,rprimd,xcart,xred)

 call var2hist(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,zDEBUG)

 hist%histV(:,:,hist%ihist)=vel(:,:)
 hist%histT(hist%ihist)=itime*ab_mover%dtion

end subroutine pred_velverlet
!!***

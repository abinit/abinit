!{\src2tex{textfont=tt}}
!!****f* ABINIT/pred_hmc
!! NAME
!!  pred_hmc
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


subroutine pred_hmc(ab_mover,hist,itime,icycle,ntime,ncycle,zDEBUG,iexit)
    
 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_abimover
 use m_abihist

 implicit none

!Arguments ------------------------------------
 type(abimover),intent(in)   :: ab_mover
 type(abihist),intent(inout) :: hist
 integer,intent(in)          :: itime
 integer,intent(in)          :: icycle
 integer,intent(in)          :: ntime
 integer,intent(in)          :: ncycle
 integer,intent(in)          :: iexit
 logical,intent(in)          :: zDEBUG

!Local variables-------------------------------

 integer  :: ii,jj                                                              ! dummy integers for loop indexes
 real(dp) :: etotal,epot,ekin,ekin_tmp                                          ! total, potential (electronic), kinetic (ionic) energies  
 real(dp) :: xcart(3,ab_mover%natom),xcart_next(3,ab_mover%natom)               ! Cartesian coordinates of all ions
 real(dp) :: xred(3,ab_mover%natom),xred_next(3,ab_mover%natom)                 ! reduced coordinates of all ions
 real(dp) :: vel(3,ab_mover%natom),vel_nexthalf(3,ab_mover%natom)               ! ionic velocities in Cartesian coordinates
 real(dp) :: fcart(3,ab_mover%natom),fred(3,ab_mover%natom)                     ! forces, Cartesian and reduced coordinates
 real(dp) :: fred_corrected(3,ab_mover%natom),fcart_corrected(3,ab_mover%natom)
 real(dp) :: factor                                                             ! factor, indicating change of time step at last iteration

 real(dp) :: acell(3)                                                           ! lattice parameters
 real(dp) :: rprimd(3,3),rprim(3,3)                                             ! lattice vectors
 real(dp),allocatable,save :: vel_prev(:,:)                                     ! velocities at the end of each time step (half time step ahead of coordinates)

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


 if(icycle==1)then
 !get current value of potential energy

 !generate new set of velocities

 !save current values of energies, as well as xred,xcart,etc.

 elseif(icycle==12)then

 !get the new values of energies

 !check if the new state is to be accepted

 !if the state is accepted, update the history and coordinates, otherwise retract

 elseif(icycle>12)then

 endif

end subroutine pred_hmc
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/testsusmat
!! NAME
!! testsusmat
!!
!! FUNCTION
!! test wether a new susceptibility matrix and/or a new dielectric matrix must be computed 
!! and return the logical result
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA,XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dielop : option for the computation of the dielectric matrix
!! dtset  : 
!! istep  : number of the current SCF cycle
!! OUTPUT
!! compute : 
!!  * if dielop >= 1 and istep == 1 => TRUE
!!  * if dielop >= 2 and istep == dtset%dielstrt => TRUE
!!  * if (dtset%iprcel >= 140 and <=170) depends on the periodicity modulo 10 of istep and iprcel
!!  * otherwise FALSE
!!
!! PARENTS
!!      prcref,prcref_PMA,vtorho
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine testsusmat(compute,dielop,dielstrt,dtset,istep)

 use m_profiling_abi

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'testsusmat'
!End of the abilint section

 implicit none

!Arguments-------------------------------
!scalars
 integer,intent(in) :: dielop,dielstrt,istep
 logical,intent(out) :: compute
 type(dataset_type),intent(in) :: dtset

!Local variables -------------------------

! *********************************************************************

 compute=.FALSE.
 if((dtset%iprcel >= 140).and.(dtset%iprcel<=170)) then
   if(modulo(dtset%iprcel,10).ne.0) then
     compute=(modulo(istep,modulo(dtset%iprcel,10))==0)
   else
     compute=(modulo(istep,10)==0)
   end if
 end if
 if (istep==1 .and. dielop>=2) compute=.TRUE.
 if (istep==dielstrt .and. dielop>=1) compute=.TRUE.
!DEBUG
!if (compute) then
!write(std_err,*) 'testsusmat : TRUE'
!else
!write(std_err,*) 'testsusmat : FALSE',dielop,dielstrt,istep,dtset%iprcel,modulo(istep,10),&
!&modulo(dtset%iprcel,10),modulo(dtset%iprcel,modulo(dtset%iprcel,10))
!end if
!ENDDEBUG
end subroutine testsusmat


 
!if( (istep==1        .and. dielop>=2) .or. &
!&     (istep==dielstrt .and. dielop>=1) .or. &
!&       computesusmat       )then

!if((iprcel >= 140) .and. (iprcel <= 170)) then
!if(modulo(iprcel,10).ne.0) then
!computediel=(modulo(istep,10)==modulo(iprcel,modulo(iprcel,10))) 
!else
!computediel=(modulo(istep,10)==0)
!end if
!end if
!
!if( (istep==1        .and. dielop>=2) &
!&     .or. (istep==dielstrt .and. dielop>=1) &
!&     .or. computediel          )then
!!***

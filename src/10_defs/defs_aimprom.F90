!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_aimprom
!! NAME
!! defs_aimprom
!!
!! FUNCTION
!! Default declarations, and interfaces for the aim utility.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_aimprom

 use defs_basis

 implicit none

 ! UNITS
 integer, save :: unt0,unto,unt,untc,unts,untd,untl,untg,unta,untad,untp,untout
 integer,save :: aim_iomode
 ! DRIVER VARIABLES
 real(dp), save :: maxatdst,maxcpdst
 integer, parameter :: ndif=45,ngaus=200,npos=1000 
 integer, allocatable, save :: typat(:), corlim(:)
 integer, save :: ntypat,nnpos,natom
 integer, save :: nsimax,batcell,npc,nbcp,nrcp,nccp
 integer, save :: icpc(npos*ndif),npcm3,slc
 real(dp), save :: rprimd(3,3),ivrprim(3,3),trivrp(3,3)
 real(dp), allocatable, save :: xred(:,:),xatm(:,:),rminl(:)
 real(dp), save :: tpi,sqfp,fpi,sqpi,sqtpi,atp(3,npos)
 real(dp), save :: h0,hmin,r0,ttsrf,ttcp,tttot
 real(dp), save :: cth(ngaus),th(ngaus),ph(ngaus),wcth(ngaus),wph(ngaus),rs(ngaus,ngaus)
 real(dp), save :: pc(3,npos*ndif), evpc(3,npos*ndif),zpc(3,3,npos*ndif), pcrb(3,npos*ndif)
 logical, save :: deb,ldeb
!!! interface chgbas
!!!    subroutine bschg1(vv,dir)
!!!      implicit none
!!!      integer, intent(in) :: dir
!!!      real(dp),intent(inout) :: vv(3)
!!!    end subroutine bschg1
!!!    subroutine bschg2(aa,dir)
!!!      implicit none
!!!      integer, intent(in) :: dir
!!!      real(dp),intent(inout) :: aa(3,3)
!!!    end subroutine bschg2
!!! end interface chgbas

 contains 
!!***

!----------------------------------------------------------------------

!!****f* defs_aimprom/aim_shutdown
!! NAME
!!  aim_shutdown
!!
!! FUNCTION
!!  Free memory allocated in the module. Close units. Mainly used to pass the abirules
!!
!! PARENTS
!!      aim
!!
!! CHILDREN
!!
!! SOURCE

 subroutine aim_shutdown()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'aim_shutdown'
!End of the abilint section

  implicit none

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'aim_shutdown'
!End of the abilint section

!Local variables-------------------------------
 integer :: ii
 logical :: is_open
 integer :: all_units(12) 

 ! *********************************************************************

 !if (allocated(typat)) then
 !  ABI_FREE(typat)
 !end if
 !if (allocated(corlim)) then
 !  ABI_FREE(corlim)
 !end if
 !if (allocated(xred)) then
 !  ABI_FREE(xred)
 !end if
 !if (allocated(xatm)) then
 !  ABI_FREE(rminl)
 !end if

 all_units(:) = [unt0,unto,unt,untc,unts,untd,untl,untg,unta,untad,untp,untout]
 do ii=1,size(all_units)
   inquire(unit=all_units(ii), opened=is_open)
   if (is_open) close(all_units(ii))
 end do

 end subroutine aim_shutdown
!!***

end module defs_aimprom
!!***

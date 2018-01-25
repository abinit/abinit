!{\src2tex{textfont=tt}}
!!****f* ABINIT/xfh_update
!! NAME
!! xfh_update
!!
!! FUNCTION
!! Update the contents of the history xfhist taking values
!! from xred, acell, rprim, fred_corrected and strten
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, JCC, SE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!!
!! SIDE EFFECTS
!!
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine xfh_update(ab_xfh,acell,fred_corrected,natom,rprim,strten,xred)

 use m_profiling_abi
 use defs_basis
 use m_abimover

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xfh_update'
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars

type(ab_xfh_type),intent(inout) :: ab_xfh
integer,intent(in) :: natom

!arrays
real(dp),intent(in) :: acell(3)
real(dp),intent(in) :: xred(3,natom)
real(dp),intent(in) :: rprim(3,3)
real(dp),intent(in) :: fred_corrected(3,natom)
real(dp),intent(in) :: strten(6)

!Local variables-------------------------------
!scalars
!integer :: kk

!*********************************************************************

!DEBUG
!write (ab_out,*) '---WROTE TO XFHIST---'

!write (ab_out,*) 'XRED'
!do kk=1,natom
!write (ab_out,*) xred(:,kk)
!end do
!write (ab_out,*) 'FRED'
!do kk=1,natom
!write (ab_out,*) fred_corrected(:,kk)
!end do
!write(ab_out,*) 'RPRIM'
!do kk=1,3
!write(ab_out,*) rprim(:,kk)
!end do
!write(ab_out,*) 'ACELL'
!write(ab_out,*) acell(:)
!DEBUG

 ab_xfh%nxfh=ab_xfh%nxfh+1

 ab_xfh%xfhist(:,1:natom,1,ab_xfh%nxfh)=xred(:,:)
 ab_xfh%xfhist(:,natom+1,1,ab_xfh%nxfh)=acell(:)
 ab_xfh%xfhist(:,natom+2:natom+4,1,ab_xfh%nxfh)=rprim(:,:)
 ab_xfh%xfhist(:,1:natom,2,ab_xfh%nxfh)=fred_corrected(:,:)
 ab_xfh%xfhist(:,natom+2,2,ab_xfh%nxfh)=strten(1:3)
 ab_xfh%xfhist(:,natom+3,2,ab_xfh%nxfh)=strten(4:6)

end subroutine xfh_update
!!***

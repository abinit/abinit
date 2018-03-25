!{\src2tex{textfont=tt}}
!!****f* ABINIT/inspln
!! NAME
!! inspln
!!
!! FUNCTION
!! This procedure gives the values of the spline coefficients
!! (second derivatives) in the 1D grid with periodic boundary
!! conditions at rsid - the values of the unknown functions specified
!! in the vector valf of direction idir
!!
!! COPYRIGHT
!! Copyright (C) 2002-2018 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  idir= direction following which the derivatives are evaluated
!!  snn, tnn=remaining bi-dimensional coordinates of the line along which
!!        the derivative is to be computed
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  This routine works on the data contained in the aimfields module
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      initaim
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine inspln(idir,snn,tnn)

 use defs_basis
 use defs_aimfields
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inspln'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,snn,tnn

!Local variables-------------------------------
!scalars
 integer :: dim,ii
 real(dp) :: ss
!arrays
 real(dp) :: rsid(ngfft(idir)),valf(ngfft(idir))
 real(dp),pointer :: ptc(:),ptd(:),ptp(:)

! *************************************************************************

!POINTER INITIALIZATION

 if (idir==1) then
   valf(:)=dvl(:,snn,tnn)
 elseif (idir==2) then
   valf(:)=dvl(tnn,:,snn)
 else
   valf(:)=dvl(snn,tnn,:)
 end if

 nullify(ptd,ptc,ptp)
 if(idir==1) then
   ptd=>dig1;ptc=>cdig1;ptp=>llg1
 elseif (idir==2) then
   ptd=>dig2;ptc=>cdig2;ptp=>llg2
 else
   ptd=>dig3;ptc=>cdig3;ptp=>llg3
 end if

 dim=ngfft(idir)

!FIRST CYCLE OF RECURRENCE

 rsid(1)=valf(2)+valf(dim)-2.*valf(1)
 rsid(1)=rsid(1)/ptd(1)
 do ii=2,dim-1
   rsid(ii)=valf(ii+1)+valf(ii-1)-2.*valf(ii)
   rsid(ii)=(rsid(ii)-ptc(ii-1)*rsid(ii-1))/ptd(ii)
 end do
 ss=0._dp
 do ii=1,dim-1
   ss=ss+rsid(ii)*ptp(ii)
 end do
 rsid(dim)=valf(1)+valf(dim-1)-2.*valf(dim)
 rsid(dim)=(rsid(dim)-ss)/ptd(dim)

!SECOND CYCLE WITH TRANSPOSED MATRIX

 rsid(dim)=rsid(dim)/ptd(dim)
 rsid(dim-1)=(rsid(dim-1)-ptc(dim-1)*rsid(dim))/ptd(dim-1)
 do ii=dim-2,1,-1
   rsid(ii)=(rsid(ii)-ptc(ii)*rsid(ii+1)-ptp(ii)*rsid(dim))/ptd(ii)
 end do

 if (idir==1) then
   ddx(:,snn,tnn)=rsid(:)
 elseif (idir==2) then
   ddy(tnn,:,snn)=rsid(:)
 else
   ddz(snn,tnn,:)=rsid(:)
 end if

end subroutine inspln
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp8cc
!! NAME
!! psp8cc
!!
!! FUNCTION
!! Compute the core charge density, for use in the XC core
!! correction, following the function definition valid
!! for format 8 of the pseudopotentials.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mmax=maximum number of points in real space grid in the psp file
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  rchrg=cut-off radius for the core density
!!
!! OUTPUT
!!  xccc1d(n1xccc,6)= 1D core charge function and its four first derivatives
!!
!! PARENTS
!!      psp8in
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp8cc(mmax,n1xccc,rchrg,xccc1d)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp8cc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mmax,n1xccc
 real(dp),intent(in) :: rchrg
!arrays
 real(dp),intent(inout) :: xccc1d(n1xccc,6) !vz_i

!Local variables-------------------------------
!scalars
 integer :: i1xccc,idum,irad,jj
 real(dp) :: amesh,c1,c2,c3,c4,damesh,dri,pi4i,tff,xp,xpm1,xpm2,xpp1,xx
 character(len=500) :: message,errmsg
!arrays
 real(dp) :: rscale(5)
 real(dp),allocatable :: ff(:,:),rad(:)

!**********************************************************************

 ABI_ALLOCATE(ff,(mmax,5))
 ABI_ALLOCATE(rad,(mmax))

 pi4i=quarter/pi
!
!Read from pp file the model core charge and its first 4 derivatives
!assumed to be on a linear grid starting at zero.
!The input functions contain the 4pi factor, and must be rescaled.

 do irad=1,mmax
   read(tmp_unit,*, err=10, iomsg=errmsg) idum,rad(irad),(ff(irad,jj),jj=1,5)
 end do


!Check that rad grid is linear starting at zero
 amesh=rad(2)-rad(1)
 damesh=zero
 do irad=2,mmax-1
   damesh=max(damesh,abs(rad(irad)+amesh-rad(irad+1)))
 end do

 if(damesh>tol8 .or. rad(1)/=zero) then
   write(message, '(5a)' )&
&   'Pseudopotential input file requires linear radial mesh',ch10,&
&   'starting at zero.',ch10,&
&   'Action : check your pseudopotential input file.'
   MSG_ERROR(message)
 end if

!Check that input rchrg is consistent with last grid point
 if(rchrg>rad(mmax)) then
   write(message, '(5a)' )&
&   'Pseudopotential input file core charge mesh',ch10,&
&   'is inconsistent with rchrg in header.',ch10,&
&   'Action : check your pseudopotential input file.'
   MSG_ERROR(message)
 end if

!Factors for unit range scaling
 do jj = 1, 5
   rscale(jj)=rchrg**(jj-1)
 end do

!Generate uniform mesh xx in the box cut by rchrg
!and interpolate the core charge and derivatives
!Cubic polynomial interpolation is used which is consistent
!with the original interpolation of these functions from
!a log grid to the input linear grid.

 dri=1.d0/amesh
 do i1xccc=1,n1xccc
   xx=(i1xccc-1)* rchrg/dble(n1xccc-1)

!  index to find bracketing input mesh points
   irad = int(dri * xx) + 1
   irad = max(irad,2)
   irad = min(irad,mmax-2)
!  interpolation coefficients
   xp = dri * (xx - rad(irad))
   xpp1 = xp + one
   xpm1 = xp - one
   xpm2 = xp - two
   c1 = -xp * xpm1 * xpm2 * sixth
   c2 = xpp1 * xpm1 * xpm2 * half
   c3 = - xp * xpp1 * xpm2 * half
   c4 = xp * xpp1 * xpm1 * sixth
!  Now do the interpolation on all derivatives for this grid point
!  Include 1/4pi normalization and unit range scaling
   do jj=1,5
     tff =  c1 * ff(irad - 1, jj) &
&     + c2 * ff(irad    , jj) &
&     + c3 * ff(irad + 1, jj) &
&     + c4 * ff(irad + 2, jj)
     xccc1d(i1xccc,jj)=pi4i*rscale(jj)*tff
   end do
 end do

!5th derivative is apparently not in use, so set to zero
 xccc1d(:,6)=zero

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(rad)

 return 

 ! Handle IO error
 10 continue
 MSG_ERROR(errmsg)

end subroutine psp8cc
!!***

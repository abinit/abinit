!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp9cc
!! NAME
!! psp9cc
!!
!! FUNCTION
!! Compute the core charge density, for use in the XC core
!! correction, following the function definition valid
!! for format 9 of the pseudopotentials (PSML).
!!
!! COPYRIGHT
!! Copyright (C) 2017 ABINIT group (YP)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mmax=maximum number of points in real space grid in the psp file
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!
!! OUTPUT
!!  rchrg=cut-off radius for the core density
!!  xccc1d(n1xccc,6)= 1D core charge function and its four first derivatives
!!
!! NOTES
!!  This routine will be built only if PSML support is enabled.
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

subroutine psp9cc(psxml,mmax,n1xccc,rad,rchrg,xccc1d)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_psml

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp9cc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mmax,n1xccc
 real(dp),intent(out) :: rchrg
 type(ps_t),intent(in) :: psxml
!arrays
 real(dp),intent(in) :: rad(mmax)
 real(dp),intent(inout) :: xccc1d(n1xccc,6) !vz_i

!Local variables-------------------------------
!scalars
 integer :: i1xccc,idum,irad,jj
 real(dp) :: amesh,c1,c2,c3,c4,damesh,dri,pi4i,tff,xp,xpm1,xpm2,xpp1,xx,twelvth
 character(len=500) :: message,errmsg
!arrays
 integer :: iwork(8)
 real(dp) :: rscale(5),dpoly(6,6),vpoly(6)
 real(dp),allocatable :: ff(:,:)

!**********************************************************************

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

 ABI_ALLOCATE(ff,(mmax,5))

 dri = one / amesh
 pi4i = quarter / pi
 twelvth = one / 12.0_dp

!Read from pp file the model core charge and calculate its first 4 derivatives
!assumed to be on a linear grid starting at zero.
!The input functions contain the 4pi factor, and must be rescaled.

!Store the value of the pseudo-core charge.
 ff(:,:) = zero
 do jj=1,mmax
   ff(jj,1) = ps_CoreCharge_Value(psxml,rad(jj))
 end do

!Calculate 4 first derivatives with 5-point stencil, except borders
 do irad=3,mmax-2
   ff(irad,2) = (-ff(irad+2,1) + 8.0d0*ff(irad+1,1) - &
&    8.0d0*ff(irad-1,1) + ff(irad-2,1)) * twelvth * dri
   ff(irad,3) = (-ff(irad+2,1) + 16.0d0*ff(irad+1,1) - 30.0d0*ff(irad,1) + &
&    16.0d0*ff(irad-1,1) - ff(irad-2,1)) * twelvth * dri * dri
   ff(irad,4) = (ff(irad+2,1) - 2.0d0*ff(irad+1,1) + &
&    2.0d0*ff(irad-1,1) - ff(irad-2,1)) * half * dri * dri * dri
   ff(irad,5) = (ff(irad+2,1) - 4.0d0*ff(irad+1,1) + 6.0d0*ff(irad,1) - &
&    4.0d0*ff(irad-1,1) + ff(irad-2,1)) * dri * dri * dri * dri
 end do

!Add border near zero using polynomial fit
 dpoly(:,:) = zero
 dpoly(:,1) = one
 vpoly(:) = zero
 vpoly(1) = ff(1,1)
 do irad=2,6
   do jj=1,6
     dpoly(irad,jj) = rad(irad)**(jj-1)
   end do
   vpoly(irad) = ff(irad,1)
 end do
 call dgesv(6,1,dpoly,6,iwork,vpoly,6,idum)

 do irad=1,2
   ff(irad,2) = &
&    vpoly(2) + 2.0d0*vpoly(3)*rad(irad) + &
&    3.0d0*vpoly(4)*rad(irad)*rad(irad) + &
&    4.0d0*vpoly(5)*rad(irad)*rad(irad)*rad(irad) + &
&    5.0d0*vpoly(6)*rad(irad)*rad(irad)*rad(irad)*rad(irad)
   ff(irad,3) = &
&    2.0d0*vpoly(3)*rad(irad) + &
&    6.0d0*vpoly(4)*rad(irad) + &
&    12.0d0*vpoly(5)*rad(irad)*rad(irad) + &
&    20.0d0*vpoly(6)*rad(irad)*rad(irad)*rad(irad)
   ff(irad,4) = &
&    6.0d0*vpoly(4) + &
&    24.0d0*vpoly(5)*rad(irad) + &
&    60.0d0*vpoly(6)*rad(irad)*rad(irad)
   ff(irad,5) = 24.0d0*vpoly(5) + &
&    120.0d0*vpoly(6)*rad(irad)
 end do

!Make linear approximation for the tail near mmax
 do irad=1,2
   ff(mmax-2+irad,2) = ff(mmax-2,2) + irad * (ff(mmax-2,2) - ff(mmax-3,2))
   ff(mmax-2+irad,3) = ff(mmax-2,3) + irad * (ff(mmax-2,3) - ff(mmax-3,3))
   ff(mmax-2+irad,4) = ff(mmax-2,4) + irad * (ff(mmax-2,4) - ff(mmax-3,4))
   ff(mmax-2+irad,5) = ff(mmax-2,5) + irad * (ff(mmax-2,5) - ff(mmax-3,5))
 end do

!Renormalize core charge
 ff(:,:) = ff(:,:) * pi4i

!determine xcccrc where the pseudocore becomes 0
!This is a difference with respect the Hamann's treatment of the core
!charge when reading PSP8.
!In Hamann's case (PSP8), xcccrc = rchrg, and this value is
!introduced in the pseudopotential input file.
!rchrg is not included in the PSML format
 rchrg = zero
 do jj=mmax,1,-1
   if (ff(jj,1) > tol13) then
     rchrg=rad(jj)
     exit
   end if
 end do

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

end subroutine psp9cc
!!***

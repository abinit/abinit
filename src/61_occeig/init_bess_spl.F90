!{\src2tex{textfont=tt}}
!!****f* ABINIT/init_bess_spl
!! NAME
!! init_bess_spl
!!
!! FUNCTION
!! Pre-calculate the j_v(y) for recip_ylm on regular grid
!!     NOTE: spherical Bessel function small j!
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  mbess=  max number of points on grid for integral
!!  (THIS VARIABLE WAS UNUSED AND REMOVED DURING PEAUTIFICATION. PMA. ) bessargmax= max point to which we will integrate
!!  bessint_delta = space between integral arguments
!!  mlang=  max angular momentum
!!
!! OUTPUT
!!  bess_spl=array of integrals
!!  bess_spl_der=array of derivatives of integrals
!!  x_bess=coordinates of points belonging to the grid
!!
!! PARENTS
!!      m_cut3d,partial_dos_fractions
!!
!! CHILDREN
!!      besjm,spline
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine init_bess_spl(mbess,bessint_delta,mlang,bess_spl,bess_spl_der,x_bess)


 use defs_basis
 use m_splines
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_bess_spl'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mbess,mlang
 real(dp),intent(in) :: bessint_delta
!arrays
 real(dp),intent(out) :: bess_spl(mbess,mlang),bess_spl_der(mbess,mlang)
 real(dp),intent(out) :: x_bess(mbess)

!Local variables -------------------------
! function calls (NumRec)
!scalars
 integer :: iintarg,ll
 character(len=500) :: msg
 real(dp) :: cosdelta_bess,sindelta_bess,yp1,ypn
!arrays
 real(dp),allocatable :: cosbessx(:),sinbessx(:)

! *********************************************************************

!DEBUG
!write(std_out,*) 'init_bess_spl enter '
!ENDDEBUG
 if (mbess < 2) then
   msg = 'init_bess_spl  Error: need more than one point for the interpolation routines'
   MSG_ERROR(msg)
 end if

!-----------------------------------------------------------------
!Bessel function into array
!-----------------------------------------------------------------

!integration grid is nfiner times finer than the interpolation grid
 ABI_ALLOCATE(sinbessx,(mbess))
 ABI_ALLOCATE(cosbessx,(mbess))

 sindelta_bess = sin(bessint_delta)
 cosdelta_bess = cos(bessint_delta)
!
!could be done by chain rule for cos sin (is it worth it?) but
!precision problems as numerical errors are propagated.
!
 do iintarg=1,mbess
   x_bess(iintarg) = (iintarg-1)*bessint_delta
   sinbessx(iintarg) = sin(x_bess(iintarg))
   cosbessx(iintarg) = cos(x_bess(iintarg))

!  x_bess(iintarg) = x_bess(iintarg-1)+bessint_delta
!  !  get sin and cos of x_bess arguments
!  sinbessx(iintarg) = sinbessx(iintarg-1)*cosdelta_bess &
!  &                     + cosbessx(iintarg-1)*sindelta_bess
!  cosbessx(iintarg) = cosbessx(iintarg-1)*cosdelta_bess &
!  &                     - sinbessx(iintarg-1)*sindelta_bess
 end do

!write(std_out,*) 'x_bess = ', x_bess
!
!fill bess_spl array
!
 do ll=0,mlang-1

   call besjm(one,bess_spl(:,ll+1),cosbessx,   &
&   ll,mbess,sinbessx,x_bess)
!  
!  call spline to get 2nd derivative (reuse in splint later)
!  
   yp1 = zero
   ypn = zero
   call spline (x_bess, bess_spl(:,ll+1), mbess, yp1, ypn, &
&   bess_spl_der(:,ll+1))
 end do

!DEBUG
!write(std_out,*) ' bess funct  0   1   2   3   4'
!do iintarg=1,mbess
!write(std_out,*) x_bess(iintarg), (bess_spl(iintarg,ll),ll=1,mlang)
!end do
!ENDDEBUG

 ABI_DEALLOCATE(sinbessx)
 ABI_DEALLOCATE(cosbessx)

end subroutine init_bess_spl
!!***

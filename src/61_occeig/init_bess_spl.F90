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
!!  nx = max number of points on grid for integral
!!  (THIS VARIABLE WAS UNUSED AND REMOVED DURING PEAUTIFICATION. PMA. ) bessargmax= max point to which we will integrate
!!  delta = space between integral arguments
!!  mlang= max angular momentum
!!
!! OUTPUT
!!  bess_spl=array of integrals
!!  bess_spl_der=array of derivatives of integrals
!!  xx=coordinates of points belonging to the grid
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

subroutine init_bess_spl(nx,delta,mlang,bess_spl,bess_spl_der,xx)

 use defs_basis
 use m_splines
 use m_profiling_abi
 use m_errors

 use m_special_funcs,  only : besjm

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_bess_spl'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,mlang
 real(dp),intent(in) :: delta
!arrays
 real(dp),intent(out) :: bess_spl(nx,mlang),bess_spl_der(nx,mlang)
 real(dp),intent(out) :: xx(nx)

!Local variables -------------------------
!scalars
 integer :: ix,ll
 character(len=500) :: msg
 real(dp) :: yp1,ypn
!arrays
 real(dp),allocatable :: cosbessx(:),sinbessx(:)

! *********************************************************************

 if (nx < 2) then
   MSG_ERROR('need more than one point for the interpolation routines')
 end if

 !-----------------------------------------------------------------
 !Bessel function into array
 !-----------------------------------------------------------------
 ! integration grid is nfiner times finer than the interpolation grid
 ABI_ALLOCATE(sinbessx,(nx))
 ABI_ALLOCATE(cosbessx,(nx))

 ! could be done by chain rule for cos sin (is it worth it?) but
 ! precision problems as numerical errors are propagated.
 do ix=1,nx
   xx(ix) = (ix-1) * delta
   sinbessx(ix) = sin(xx(ix))
   cosbessx(ix) = cos(xx(ix))
 end do

 ! fill bess_spl array
 do ll=0,mlang-1
   call besjm(one,bess_spl(:,ll+1),cosbessx,ll,nx,sinbessx,xx)

   ! call spline to get 2nd derivative (reuse in splint later)
   yp1 = zero; ypn = zero
   call spline(xx, bess_spl(:,ll+1), nx, yp1, ypn, bess_spl_der(:,ll+1))
 end do

!write(std_out,*) ' bess funct  0   1   2   3   4'
!do ix=1,nx
!write(std_out,*) xx(ix), (bess_spl(ix,ll),ll=1,mlang)
!end do

 ABI_DEALLOCATE(sinbessx)
 ABI_DEALLOCATE(cosbessx)

end subroutine init_bess_spl
!!***

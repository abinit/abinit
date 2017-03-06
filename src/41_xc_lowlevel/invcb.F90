!{\src2tex{textfont=tt}}
!!****f* ABINIT/invcb
!! NAME
!! invcb
!!
!! FUNCTION
!! Compute a set of inverse cubic roots as fast as possible :
!! rspts(:)=rhoarr(:)$^\frac{-1}{3}$
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npts=number of real space points on which density is provided
!!  rhoarr(npts)=input data
!!
!! OUTPUT
!!  rspts(npts)=inverse cubic root of rhoarr
!!
!! PARENTS
!!      drivexc,gammapositron,xchcth,xcpbe,xcpositron,xctfw
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine invcb(rhoarr,rspts,npts)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'invcb'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npts
!arrays
 real(dp),intent(in) :: rhoarr(npts)
 real(dp),intent(out) :: rspts(npts)

!Local variables-------------------------------
!scalars
 integer :: ii,ipts
 real(dp),parameter :: c2_27=2.0e0_dp/27.0e0_dp,c5_9=5.0e0_dp/9.0e0_dp
 real(dp),parameter :: c8_9=8.0e0_dp/9.0e0_dp,m1thrd=-third
 real(dp) :: del,prod,rho,rhom1,rhomtrd
 logical :: test
!character(len=500) :: message

! *************************************************************************

!Loop over points : here, brute force algorithm
!do ipts=1,npts
!rspts(ipts)=sign( (abs(rhoarr(ipts)))**m1thrd,rhoarr(ipts))
!end do
!

!write(std_out,*)' invcb : rhoarr, rspts'

 rhomtrd=sign( (abs(rhoarr(1)))**m1thrd, rhoarr(1) )
 rhom1=one/rhoarr(1)
 rspts(1)=rhomtrd
 do ipts=2,npts
!  write(std_out,*)
!  write(std_out,*)rhoarr(ipts),rspts(ipts)
   rho=rhoarr(ipts)
   prod=rho*rhom1
!  If the previous point is too far ...
   if(prod < 0.01_dp .or. prod > 10._dp )then
     rhomtrd=sign( (abs(rho))**m1thrd , rho )
     rhom1=one/rho
   else
     del=prod-one
     do ii=1,5
!      Choose one of the two next lines, the last one is more accurate
!      rhomtrd=((one+third*del)/(one+two_thirds*del))*rhomtrd
       rhomtrd=((one+c5_9*del)/(one+del*(c8_9+c2_27*del)))*rhomtrd
       rhom1=rhomtrd*rhomtrd*rhomtrd
       del=rho*rhom1-one
!      write(std_out,*)rhomtrd,del
       test = del*del < 1.0e-24_dp
       if(test) exit
     end do
     if( .not. test) then
       rhomtrd=sign( (abs(rho))**m1thrd , rho )
     end if
   end if
   rspts(ipts)=rhomtrd
 end do

end subroutine invcb
!!***


!{\src2tex{textfont=tt}}
!!****f* ABINIT/findminscf
!!
!! NAME
!! findminscf
!!
!! FUNCTION
!! Compute the minimum of a function whose value
!! and derivative are known at two points,
!! using different algorithms.
!! Also deduce different quantities at this predicted
!! point, and at the two other points
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!

!! INPUTS
!! choice=1,uses a linear interpolation of the derivatives
!!       =2,uses a quadratic interpolation based on the
!!        values of the function, and the second derivative at mid-point
!! etotal_1=first value of the function
!! etotal_2=second value of the function
!! dedv_1=first value of the derivative
!! dedv_2=second value of the derivative
!! lambda_1=first value of the argument
!! lambda_2=second value of the argument
!!
!! OUTPUT
!! dedv_predict=predicted value of the derivative (usually zero,
!!  except if choice=4, if it happens that a minimum cannot be located,
!!  and a trial step is taken)
!! d2edv2_predict=predicted value of the second derivative (not if choice=4)
!! d2edv2_1=first value of the second derivative (not if choice=4)
!! d2edv2_2=second value of the second derivative (not if choice=4)
!! etotal_predict=predicted value of the function
!! lambda_predict=predicted value of the argument
!! status= 0 if everything went normally ;
!!         1 if negative second derivative
!!         2 if some other problem
!!
!! PARENTS
!!      scfcge
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine findminscf(choice,dedv_1,dedv_2,dedv_predict,&
& d2edv2_1,d2edv2_2,d2edv2_predict,&
& etotal_1,etotal_2,etotal_predict,&
& lambda_1,lambda_2,lambda_predict,errid,errmess)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'findminscf'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice
 integer,intent(out) :: errid
 character(len=500), intent(out) :: errmess
 real(dp),intent(in) :: dedv_1,dedv_2,etotal_1,etotal_2,lambda_1,lambda_2
 real(dp),intent(out) :: d2edv2_1,d2edv2_2,d2edv2_predict,dedv_predict
 real(dp),intent(out) :: etotal_predict,lambda_predict

!Local variables-------------------------------
!scalars
 real(dp) :: cc,d2edv2_mid,d_lambda,dedv_2bis
 real(dp) :: dedv_mid2,etotal_2bis
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(std_out,*)' findmin : enter'
!write(std_out,*)' choice,lambda_1,lambda_2=',choice,lambda_1,lambda_2
!ENDDEBUG

 errid = AB7_NO_ERROR
 d_lambda=lambda_1-lambda_2

 if(choice==1) then

!  Use the derivative information to predict lambda
   d2edv2_mid=(dedv_1-dedv_2)/d_lambda
   lambda_predict=lambda_2-dedv_2/d2edv2_mid
   dedv_predict=dedv_2+(lambda_predict-lambda_2)*d2edv2_mid
   d2edv2_1=d2edv2_mid
   d2edv2_2=d2edv2_mid
   d2edv2_predict=d2edv2_mid
!  also use the first energy to predict new energy
   etotal_predict=etotal_1+dedv_1*(lambda_predict-lambda_1)&
&   +0.5_dp*d2edv2_1*(lambda_predict-lambda_1)**2
   etotal_2bis=etotal_1+dedv_1*(lambda_2-lambda_1)&
&   +0.5_dp*d2edv2_1*(lambda_2-lambda_1)**2

   if(d2edv2_mid<0.0_dp)then
     errid = AB7_ERROR_MIXING_INTERNAL
     write(errmess,'(a,es18.10,a)')'The second derivative is negative, equal to ',d2edv2_mid,'.'
     MSG_WARNING(errmess)
   end if

 else if(choice==2) then

!  Use energies and first derivative information
!  etotal = aa + bb * lambda + cc * lambda**2
   dedv_mid2=(etotal_1-etotal_2)/d_lambda
   cc=(dedv_1-dedv_mid2)/d_lambda
   lambda_predict=lambda_1-0.5_dp*dedv_1/cc
   d2edv2_1=2*cc
   d2edv2_2=d2edv2_1
   d2edv2_predict=d2edv2_1
   if(d2edv2_predict<0.0_dp)then
     errid = AB7_ERROR_MIXING_INTERNAL
     write(errmess, '(a,es18.10,a,a,a)' )&
&     'The second derivative is negative, equal to',d2edv2_predict,'.',ch10,&
&     '=> Pivoting                     '
     MSG_WARNING(errmess)
     if(etotal_2 < etotal_1)then
       lambda_predict=lambda_2-0.5_dp*(lambda_1-lambda_2)
     else
       lambda_predict=lambda_1-0.5_dp*(lambda_2-lambda_1)
     end if
   end if
   dedv_predict=dedv_1+(lambda_predict-lambda_1)*d2edv2_1
   dedv_2bis=dedv_1+(lambda_2-lambda_1)*d2edv2_1
   etotal_predict=etotal_1+dedv_1*(lambda_predict-lambda_1)&
&   +0.5_dp*d2edv2_1*(lambda_predict-lambda_1)**2

 end if

 write(message, '(a,es12.4,a,es18.10)' ) &
& ' findmin : lambda_predict ',lambda_predict,' etotal_predict ',etotal_predict
 call wrtout(std_out,message,'COLL')

end subroutine findminscf
!!***

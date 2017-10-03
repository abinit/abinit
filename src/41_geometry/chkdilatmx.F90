!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkdilatmx
!! NAME
!! chkdilatmx
!!
!! FUNCTION
!! Check whether the new rprimd does not give a too large number
!! of plane waves, compared to the one booked for rprimd, taking
!! into account the maximal dilatation dilatmx. Actually check whether
!! the new Fermi sphere is inside the old one, dilated.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  chkdilatmx_ = if 1, will prevent to have any vector outside the Fermi sphere, possibly
!!       by rescaling (three times at most), and then stopping the execution
!!                if 0, simply send a warning, but continues execution
!!  dilatmx     = maximal dilatation factor (usually the input variable)
!!  rprimd      = new primitive vectors
!!  rprimd_orig = original primitive vectors (usually the input variable)
!!
!! OUTPUT
!!  dilatmx_errmsg=Emptry string if calculation can continue.
!!            If the calculation cannot continue, dilatmx_errmsg will contain
!!            the message that should be reported in the output file.
!!
!!            Client code should handle a possible problem with the following test:
!!
!!              if (LEN_TRIM(dilatmx_errmsg) then 
!!                dump dilatmx_errmsg to the main output file.
!!                handle_error
!!              end if 
!!    
!!
!! PARENTS
!!      driver,mover
!!
!! CHILDREN
!!      matr3eigval,matr3inv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine chkdilatmx(chkdilatmx_,dilatmx,rprimd,rprimd_orig,dilatmx_errmsg)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkdilatmx'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chkdilatmx_
 real(dp),intent(in) :: dilatmx
 character(len=500),intent(out) :: dilatmx_errmsg
!arrays
 real(dp),intent(inout) :: rprimd(3,3)
 real(dp),intent(in) :: rprimd_orig(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,mu
 real(dp) :: dilatmx_new
!arrays
 real(dp) :: eigval(3),gprimd_orig(3,3),met(3,3),old_to_new(3,3)
 real(dp) :: eigval_orig(3), alpha

! *************************************************************************

!Generates gprimd
 call matr3inv(rprimd_orig,gprimd_orig)

!Find the matrix that transform an original xcart to xred, then to the new xcart
 do mu=1,3
   old_to_new(mu,:)=rprimd(mu,1)*gprimd_orig(:,1)+&
&   rprimd(mu,2)*gprimd_orig(:,2)+&
&   rprimd(mu,3)*gprimd_orig(:,3)
 end do

!The largest increase in length will be obtained thanks
!to the diagonalization of the corresponding metric matrix :
!it is the square root of its largest eigenvalue.
 do ii=1,3
   do jj=1,3
     met(ii,jj)=old_to_new(1,ii)*old_to_new(1,jj)+&
&     old_to_new(2,ii)*old_to_new(2,jj)+&
&     old_to_new(3,ii)*old_to_new(3,jj)
   end do
 end do

 call matr3eigval(eigval,met)

 dilatmx_new=sqrt(maxval(eigval(:)))

 dilatmx_errmsg = ""
 if(dilatmx_new>dilatmx+tol6)then

! MJV 2014 07 22: correct rprim to maximum jump allowed by dilatmx
! eigenvalues of old metric tensor is needed
   do mu=1,3
     old_to_new(mu,:)=rprimd_orig(mu,1)*gprimd_orig(:,1)+&
&     rprimd_orig(mu,2)*gprimd_orig(:,2)+&
&     rprimd_orig(mu,3)*gprimd_orig(:,3)
   end do

   do ii=1,3
     do jj=1,3
       met(ii,jj)=old_to_new(1,ii)*old_to_new(1,jj)+&
&       old_to_new(2,ii)*old_to_new(2,jj)+&
&       old_to_new(3,ii)*old_to_new(3,jj)
     end do
   end do
   call matr3eigval(eigval_orig,met)
   dilatmx_new=sqrt(maxval(eigval_orig(:)))

   if(chkdilatmx_/=0)then
     alpha = (dilatmx**2 - maxval(eigval)) / (maxval(eigval_orig) - maxval(eigval))
!    for safety, only 90 percent of max jump
     alpha = 0.9_dp * alpha

     rprimd = alpha * rprimd + (one - alpha) * rprimd_orig

     write(dilatmx_errmsg,'(3a,f5.2,4a,f6.2,a)')&
&     'The new primitive vectors rprimd (an evolving quantity)',ch10,&
&     'are too large with respect to the old rprimd and the accompanying dilatmx:',dilatmx,ch10,&
&     'This large change of unit cell parameters is not allowed by the present value of dilatmx.',ch10,&
&     'Calculation continues with limited jump of ', 100._dp * alpha,' percent of the projected move.'
   else
     write(message, '(3a,es16.6,2a,es16.6,2a)' )&
&     'The new primitive vectors rprimd (an evolving quantity)',ch10,&
&     'are too large, given the initial rprimd and the accompanying dilatmx:',dilatmx,ch10,&
&     'An adequate value would have been dilatmx_new=',dilatmx_new,ch10,&
&     'As chkdilatmx=1, assume experienced user. Execution will continue.'
     MSG_WARNING(message)

 end if

end subroutine chkdilatmx
!!**

!{\src2tex{textfont=tt}}
!!****f* ABINIT/erlxconv
!! NAME
!!  erlxconv
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2016-2017 ABINIT group (DW)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine erlxconv(etotals,iexit,ihist,itime,mxhist,ntime,tolmxde)
    
 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_fstrings,  only : indent

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'erlxconv'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ihist,itime,mxhist,ntime
 integer,intent(inout) :: iexit
 real(dp), intent(in) :: tolmxde
!arrays
 real(dp),intent(in) :: etotals(mxhist)

!Local variables-------------------------------
 real(dp) :: ediff1,ediff2,maxediff
 character(len=500) :: message                   
 
! *************************************************************************

 if (ihist<3) then
   write(message, '(a,a,a)' ) ch10,&
&   ' erlxconv : minimum 3 Broyd/MD steps to check convergence of energy in relaxations',ch10
   call wrtout(std_out,message,'COLL')
 else
   ediff1 = etotals(ihist)-etotals(ihist-1)
   ediff2 = etotals(ihist)-etotals(ihist-2)
   if ((abs(ediff1)<tolmxde).and.(abs(ediff2)<tolmxde)) then
     write(message, '(a,a,i4,a,a,a,a,a,es11.4,a,a)' ) ch10,&
&     ' At Broyd/MD step',itime,', energy is converged : ',ch10,&
&     '  the difference in energy with respect to the two ',ch10,&
&     '  previous steps is < tolmxde=',tolmxde,' ha',ch10
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     iexit=1
   else
     maxediff = max(abs(ediff1),abs(ediff2))
     if(iexit==1)then
       write(message, '(a,a,a,a,i5,a,a,a,es11.4,a,es11.4,a,a)' ) ch10,&
&       ' erlxconv : WARNING -',ch10,&
&       '  ntime=',ntime,' was not enough Broyd/MD steps to converge energy: ',ch10,&
&       '  max difference in energy =',maxediff,' > tolmxde=',tolmxde,' ha',ch10
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')

       write(std_out,"(8a)")ch10,&
&       "--- !RelaxConvergenceWarning",ch10,&
&       "message: | ",ch10,TRIM(indent(message)),ch10,&
&       "..."
     else
       write(message, '(a,a,i4,a,a,a,es11.4,a,es11.4,a,a)' ) ch10,&
&       ' erlxconv : at Broyd/MD step',itime,', energy has not converged yet. ',ch10,&
&       '  max difference in energy=',maxediff,' > tolmxde=',tolmxde,' ha',ch10
       call wrtout(std_out,message,'COLL')
     end if
   end if
 end if

end subroutine erlxconv
!!***

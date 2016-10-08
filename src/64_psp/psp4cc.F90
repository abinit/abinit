!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp4cc
!! NAME
!! psp4cc
!!
!! FUNCTION
!! Compute the core charge density, for use in the XC core
!! correction, following the function definition valid
!! for the format 4 of pseudopotentials.
!! This is a even polynomial of 24th order for core density,
!! that is cut off exactly beyond rchrg.
!! It has been produced on 7 May 1992 by M. Teter.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (XG, DCA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  fchrg=magnitude of the core charge correction
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!
!! OUTPUT
!!  xccc1d(n1xccc,6)= 1D core charge function and its five first derivatives
!!
!! NOTES
!! The argument of xccc1d is assumed to be normalized, and to vary
!! from xx=0 to 1 (from r=0 to r=xcccrc)
!!
!! WARNINGS
!! the fifth derivative is not yet delivered.
!!
!! PARENTS
!!      psp1in
!!
!! CHILDREN
!!      spline
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp4cc(fchrg,n1xccc,xccc1d)

 use defs_basis
 use m_profiling_abi
 use m_splines
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp4cc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1xccc
 real(dp),intent(in) :: fchrg
!arrays
 real(dp),intent(inout) :: xccc1d(n1xccc,6) !vz_i

!Local variables-------------------------------
!scalars
 integer :: i1xccc,ider
 real(dp),parameter :: a10=-0.1156854803757563d5,a12=+0.2371534625455588d5
 real(dp),parameter :: a14=-0.3138755797827918d5,a16=+0.2582842713241039d5
 real(dp),parameter :: a18=-0.1200356429115204d5,a20=+0.2405099057118771d4
 real(dp),parameter :: a2=-0.8480751097855989d1,a4=+0.9684600878284791d2
 real(dp),parameter :: a6=-0.7490894651588015d3,a8=+0.3670890998130434d4
 real(dp) :: der1,dern,factor,gpp_1,gpp_2,gpp_3
 real(dp) :: gg,gp,gpp,gpp1,gpp2,gpp3
 character(len=500) :: message
!arrays
 real(dp),allocatable :: ff(:),ff2(:),work(:),xx(:)
 real(dp) :: x

! *************************************************************************

 ABI_ALLOCATE(ff,(n1xccc))
 ABI_ALLOCATE(ff2,(n1xccc))
 ABI_ALLOCATE(work,(n1xccc))
 ABI_ALLOCATE(xx,(n1xccc))


 if(n1xccc > 1)then
   factor=1.0d0/dble(n1xccc-1)
   do i1xccc=1,n1xccc
     xx(i1xccc)=(i1xccc-1)*factor
   end do
 else
   write(message, '(a,i0)' )'  n1xccc should larger than 1, while it is n1xccc=',n1xccc
   MSG_BUG(message)
 end if

!Initialization, to avoid some problem with some compilers
 xccc1d(1,:)=zero ; xccc1d(n1xccc,:)=zero

!Take care of each derivative separately
 do ider=0,2

   if(ider==0)then
!    Generate spline fitting for the function gg
     do i1xccc=1,n1xccc
!      ff(i1xccc)=fchrg*gg(xx(i1xccc))
       ff(i1xccc)=fchrg*gg_psp4(xx(i1xccc))
     end do
!    Complete with derivatives at end points
     der1=0.0d0
!    dern=fchrg*gp(1.0d0) 
     dern=fchrg*gp_psp4(1.0d0) 
   else if(ider==1)then
!    Generate spline fitting for the function gp
     do i1xccc=1,n1xccc
!      ff(i1xccc)=fchrg*gp(xx(i1xccc))
       ff(i1xccc)=fchrg*gp_psp4(xx(i1xccc))
     end do
!    Complete with derivatives at end points, already estimated
     der1=xccc1d(1,ider+2)
     dern=xccc1d(n1xccc,ider+2)
   else if(ider==2)then
!    Generate spline fitting for the function gpp
!    (note : the function gpp has already been estimated, for the spline
!    fitting of the function gg, but it is replaced here by the more
!    accurate analytic derivative)
     do i1xccc=1,n1xccc
       x=xx(i1xccc)
       ff(i1xccc)=fchrg*(gpp_1_psp4(x)+gpp_2_psp4(x)+gpp_3_psp4(x))
!      ff(i1xccc)=fchrg*gpp(xx(i1xccc))
     end do
!    Complete with derivatives of end points
     der1=xccc1d(1,ider+2)
     dern=xccc1d(n1xccc,ider+2)
   end if

!  Produce second derivative numerically, for use with splines
   call spline(xx,ff,n1xccc,der1,dern,ff2)
   xccc1d(:,ider+1)=ff(:)
   xccc1d(:,ider+3)=ff2(:)

 end do

 xccc1d(:,6)=zero

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(ff2)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(xx)

!DEBUG
!write(std_out,*)' psp1cc : output of core charge density and derivatives '
!write(std_out,*)'   xx          gg           gp  '
!do i1xccc=1,n1xccc
!write(std_out,'(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,1),xccc1d(i1xccc,2)
!end do
!write(std_out,*)'   xx          gpp          gg2  '
!do i1xccc=1,n1xccc
!write(std_out,'(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,3),xccc1d(i1xccc,4)
!end do
!write(std_out,*)'   xx          gp2          gpp2  '
!do i1xccc=1,n1xccc
!write(std_out,'(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,5),xccc1d(i1xccc,6)
!end do
!write(std_out,*)' psp1cc : debug done, stop '
!stop
!ENDDEBUG

 contains
 
   function gg_psp4(x)
!Expression of 7 May 1992


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gg_psp4'
!End of the abilint section

   real(dp) :: gg_psp4
   real(dp),intent(in) :: x
   gg_psp4=(1.d0+x**2*(a2 +x**2*(a4 +x**2*(a6 +x**2*(a8 + &
&   x**2*(a10+x**2*(a12+x**2*(a14+x**2*(a16+ &
&   x**2*(a18+x**2*(a20)))))))))))          *(1.0d0-x**2)**2
 end function gg_psp4

 function gp_psp4(x)
!gp(x) is the derivative of gg(x) wrt x


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gp_psp4'
!End of the abilint section

   real(dp) :: gp_psp4
   real(dp),intent(in) :: x
   gp_psp4=2.d0*x*((a2+x**2*(2.d0*a4+x**2*(3.d0*a6+x**2*(              &
&   4.d0*a8+x**2*(5.d0*a10+x**2*(6.d0*a12+x**2*(                     &
&   7.d0*a14+x**2*(8.d0*a16+x**2*(9.d0*a18+x**2*(10.d0*a20))))))))))*&
&   (1.d0-x**2)**2                                                &
&   -2.0d0*(1.d0+x**2*(a2 +x**2*(a4 +x**2*(a6 +x**2*(a8 +            &
&   x**2*(a10+x**2*(a12+x**2*(a14+x**2*(a16+            &
&   x**2*(a18+x**2*a20))))))))))        *(1.0d0-x**2) )
 end function gp_psp4

   function gpp_1_psp4(x)
!gpp(x) is the second derivative of gg(x) wrt x


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gpp_1_psp4'
!End of the abilint section

   real(dp) :: gpp_1_psp4
   real(dp),intent(in) :: x
   gpp_1_psp4= ( 2.d0*a4+ x**2*(3.d0*2.d0*a6 +x**2*(               &
&   4.d0*3.d0*a8+ x**2*(5.d0*4.d0*a10+x**2*(               &
&   6.d0*5.d0*a12+x**2*(7.d0*6.d0*a14+x**2*(               &
&   8.d0*7.d0*a16+x**2*(9.d0*8.d0*a18+x**2*(               &
&   10.d0*9.d0*a20)                                        &
&   ))))))))*(2.d0*x*(1.d0-x**2))**2
 end function gpp_1_psp4

 function gpp_2_psp4(x)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gpp_2_psp4'
!End of the abilint section

   real(dp) :: gpp_2_psp4
   real(dp),intent(in) :: x
   gpp_2_psp4=(a2+x**2*(2.d0*a4+x**2*(3.d0*a6+x**2*(                 &
&   4.d0*a8 +x**2*(5.d0*a10+x**2*(6.d0*a12+x**2*(          &
&   7.d0*a14+x**2*(8.d0*a16+x**2*(9.d0*a18+x**2*(          &
&   10.d0*a20)                                             &
&   )))))))))*(1.d0-x**2)*2*(1.d0-9.d0*x**2)
 end function gpp_2_psp4

 function gpp_3_psp4(x)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gpp_3_psp4'
!End of the abilint section

   real(dp) :: gpp_3_psp4
   real(dp),intent(in) :: x
   gpp_3_psp4=(1.d0+x**2*(a2 +x**2*(a4 +x**2*(a6 +x**2*(a8 +         &
&   x**2*(a10+x**2*(a12+x**2*(a14+x**2*(a16+         &
&   x**2*(a18+x**2*a20                               &
&   ))))))))))*(1.0d0-3.d0*x**2)*(-4.d0)
 end function gpp_3_psp4
 
end subroutine psp4cc
!!***

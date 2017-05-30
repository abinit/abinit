!{\src2tex{textfont=tt}}
!!****f* ABINIT/smeared_delta
!! NAME
!! smeared_delta
!!
!! FUNCTION
!! This subroutine calculates the smeared delta that weights matrix elements: 
!! \delta (\epsilon_{kn}-\epsilon_{k+Q,n'})
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (PB, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!! eigen0(mband*nsppol) : Eigenvalues at point K 
!! eigenq(mband*nsppol)  : Eigenvalues at point K+Q
!! mband : maximum number of bands
!! smdelta : Variable controlling the smearinf scheme
!!
!! OUTPUT
!! smdfunc(mband,mband) : Smeared delta function weight corresponding to \delta(\epsilon_{n,k} - \epsilon_{n',k+Q})
!!
!! PARENTS
!!      eig2stern,eig2tot
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine smeared_delta(eigen0,eigenq,esmear,mband,smdelta,smdfunc)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'smeared_delta'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,smdelta
!arrays
 real(dp),intent(in) :: eigen0(mband),eigenq(mband),esmear
 real(dp),intent(out) :: smdfunc(mband,mband)

!Local variables-------------------------------
!tolerance for non degenerated levels
!scalars
 integer :: ii,jj
 real(dp) :: aa,dsqrpi,gauss,xx
 character(len=500) :: message

! *********************************************************************


!---------------------------------------------------------
!Ordinary (unique) smearing function
!---------------------------------------------------------

 if(smdelta==1)then

!  Fermi-Dirac
   do ii=1,mband
     do jj= 1,mband
       xx= ( eigen0(ii) - eigenq(jj) )/esmear
       smdfunc(ii,jj)=0.25_dp/esmear/(cosh(xx/2.0_dp))**2
     end do
   end do

 else if(smdelta==2 .or. smdelta==3)then

!  Cold smearing of Marzari, two values of the "a" parameter being possible
!  first value gives minimization of the bump
   if(smdelta==2)aa=-.5634
!  second value gives monotonic occupation function
   if(smdelta==3)aa=-.8165

   dsqrpi=1.0_dp/sqrt(pi)
   do ii=1,mband
     do jj=1,mband
       xx= ( eigen0(ii) - eigenq(jj) ) / esmear
       gauss=dsqrpi*exp(-xx**2)/esmear
       smdfunc(ii,jj)=gauss*(1.5_dp+xx*(-aa*1.5_dp+xx*(-1.0_dp+aa*xx)))
     end do
   end do

 else if(smdelta==4)then

!  First order Hermite-Gaussian of Paxton and Methfessel
   dsqrpi=1.0_dp/sqrt(pi)
   do ii=1,mband
     do jj=1,mband
       xx= ( eigen0(ii) - eigenq (jj) ) / esmear
       smdfunc(ii,jj)=dsqrpi*(1.5_dp-xx**2)*exp(-xx**2)/esmear
     end do
   end do

 else if(smdelta==5)then

!  Gaussian smearing
   dsqrpi=1.0_dp/sqrt(pi)
   do ii=1,mband
     do jj=1,mband
       xx= ( eigen0(ii) - eigenq (jj) ) / esmear
       smdfunc(ii,jj)=dsqrpi*exp(-xx**2)/esmear
     end do
   end do

 else
   write(message, '(a,i0,a)' )'  Smdelta= ',smdelta,' is not allowed in smdfunc'
   MSG_BUG(message)
 end if

end subroutine smeared_delta
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/kramerskronig
!! NAME
!! kramerskronig
!!
!! FUNCTION
!! check or apply the Kramers Kronig relation:
!!  Re \epsilon(\omega) = 1 + \frac{2}{\pi}
!!  \int_0^\infty d\omega' frac{\omega'}{\omega'^2 - \omega^2} Im \epsilon(\omega') 
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2018 ABINIT group (Valerio Olevano, Lucia Reining, Francesco Sottile, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nomega=number of real frequencies
!!  omega(nomega)= real frequencies 
!!  eps(nomega)= function on the frequency grid (both real and imaginary part)
!!   real part can be used to check whether the K-K relation is satisfied or not
!!  method=method used to perform the integration
!!   0= naive integration
!!   1=simpson rule 
!!  only_check= if /=0 the real part of eps is checked against the imaginary part,
!!                a final report in written but the array eps is not modified 
!!              if ==0 the real part of eps is overwritten using the
!!              results obtained using the Kramers-Kronig relation     
!!
!! OUTPUT
!!
!! Inspired to check_kramerskronig of the DP code 
!!
!! PARENTS
!!      linear_optics_paw
!!
!! CHILDREN
!!      simpson_int,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine kramerskronig(nomega,omega,eps,method,only_check)

 use defs_basis
 use m_errors

 use m_numeric_tools,   only : simpson_int

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kramerskronig'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: method,nomega,only_check
!arrays
 real(dp),intent(in) :: omega(nomega)
 complex,intent(inout) :: eps(nomega)

!Local variables-------------------------------
!scalars
 integer,save :: enough=0
 integer :: ii,ip
 real(dp) :: acc,domega,eav,kkdif,kkrms,ww,wwp
 character(len=500) :: msg
!arrays
 real(dp) :: e1kk(nomega),intkk(nomega),kk(nomega)

! *************************************************************************

!Check whether the frequency grid is linear or not 
 domega = (omega(nomega) - omega(1)) / (nomega-1)
 do ii=2,nomega
   if (ABS(domega-(omega(ii)-omega(ii-1))) > 0.001) then 
     if (only_check/=1) then
       msg="check cannot be performed since the frequency step is not constant"
       MSG_WARNING(msg)
       RETURN
     else 
       msg=' Cannot perform integration since frequency step is not constant'
       MSG_ERROR(msg)
     end if 
   end if
 end do

!Check whether omega(1) is small or not
 if (omega(1) > 0.1/Ha_eV) then 
   if (only_check/=1) then
     msg=' Check cannot be performed since first frequency on the grid > 0.1 eV'
     MSG_WARNING(msg)
     RETURN
   else 
     msg=' Cannot perform integration since first frequency on the grid > 0.1 eV'
     MSG_ERROR(msg)
   end if
 end if 

!If eps(nomega) is not 0 warn
 if (AIMAG(eps(nomega)) > 0.1 .and. enough<50) then
   enough=enough+1
   write(msg,'(a,f8.4,3a,f8.2,2a)')&
&   'Im epsilon for omega = ',omega(nomega)*Ha_eV,' eV',ch10,&
&   'is not yet zero, epsilon_2 = ',AIMAG(eps(nomega)),ch10,&
&   'Kramers Kronig could give wrong results'
   MSG_WARNING(msg)
   if (enough==50) then 
     write(msg,'(3a)')' sufficient number of WARNINGS-',ch10,' stop writing '
     call wrtout(std_out,msg,'COLL')
   end if 
 end if
 
 
!Perform Kramers-Kronig using naive integration
 select case (method)
 case (0)

   do ii=1,nomega
     ww = omega(ii)
     acc = 0.0_dp
     do ip=1,nomega
       if (ip == ii) CYCLE
       wwp = omega(ip)
       acc = acc + wwp/(wwp**2-ww**2) *AIMAG(eps(ip))
     end do
     e1kk(ii) = one + two/pi*domega* acc
   end do
   
!    Perform Kramers-Kronig using Simpson integration    
!    Simpson O(1/N^4), from NumRec in C p 134  NumRec in Fortran p 128
 case (1)

   kk=zero

   do ii=1,nomega
     ww=omega(ii)
     do ip=1,nomega
       if (ip == ii) CYCLE
       wwp = omega(ip)
       kk(ip) = wwp/(wwp**2-ww**2) *AIMAG(eps(ip))
     end do
     call simpson_int(nomega,domega,kk,intkk)
     e1kk(ii) = one + two/pi * intkk(nomega)
   end do

 case default
   write(msg,'(a,i0)')' Wrong value for method ',method
   MSG_BUG(msg)
 end select

!at this point real part is in e1kk, need to put it into eps
 do ii=1,nomega
   eps(ii)=CMPLX(e1kk(ii),AIMAG(eps(ii)))
 end do 

!Verify Kramers-Kronig
 eav   = zero
 kkdif = zero
 kkrms = zero

 do ii=1,nomega
   kkdif = kkdif + ABS(REAL(eps(ii)) - e1kk(ii))
   kkrms = kkrms + (REAL(eps(ii)) - e1kk(ii))*(REAL(eps(ii)) - e1kk(ii))
   eav = eav + ABS(REAL(eps(ii)))
 end do

 eav = eav/nomega
 kkdif = (kkdif/nomega) / eav
 kkrms = (kkrms/nomega) / (eav*eav)

 kk = ABS(REAL(eps(1)) - e1kk(1)) / REAL(eps(1))

!Write data
 write(msg,'(a,f7.2,a)')' Kramers-Kronig transform is verified within ',MAXVAL(kk)*100,"%"
 call wrtout(std_out,msg,'COLL')

end subroutine kramerskronig
!!***

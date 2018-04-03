!{\src2tex{textfont=tt}}
!!****f* ABINIT/fsumrule
!! NAME
!! fsumrule
!!
!! FUNCTION
!! This routine performs a check of the f-sum rule:
!!
!!  \int_0^\infty d\omega \omega Im \epsilon_GGp(q,\omega) =
!!       = \frac{1}{2} \pi \omega_p^2  \frac{\rho(G-Gp)}{\rho(0)} 
!!         versor(q+G) \dot versor(q+Gp)
!!  
!!  for q = G = Gp = 0, it reads:
!!  \int_0^\infty d\omega \omega Im \epsilon_00(q=0,\omega) = \pi \omega_p^2 / 2
!!
!!  check only the second one:
!!  calculate the integral to evaluate an omega_plasma^eff to compare with omega_plasma
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2018 ABINIT group (MG,VO,LR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nomega=number of real frequencies
!!  omega(nomega)= real frequencies 
!!  eps(nomega)= function on the frequency grid 
!!   in the present implementation must be the imaginary part of epsilon_00 (q=0)
!!  method=method used to perform the integration
!!   0= naive integration
!!   1=simpson rule 
!!
!! OUTPUT
!!  Only check.
!!
!! NOTES
!!  Inspired to a similar routine of the dp code.
!!
!! PARENTS
!!
!! CHILDREN
!!      simspon_int,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fsumrule(nomega,omega,eps,omegaplasma,method)

 use m_profiling_abi

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fsumrule'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: method,nomega
 real(dp),intent(in) :: omegaplasma
!arrays
 real(dp),intent(in) :: eps(nomega),omega(nomega)

!Local variables-------------------------------
!scalars
 integer :: ii,ip
 real(dp) :: acc,check,domega,omegaplasmaeff,ww,wwp
 character(len=500) :: msg
!arrays
 real(dp) :: fs(nomega),fsint(nomega)

! *************************************************************************
 
!Check whether the frequency grid is linear or not 
 domega = (omega(nomega) - omega(1))/(nomega-1)
 do ii =2,nomega
  if (ABS(domega-(omega(ii)-omega(ii-1))) > 0.001) then
   msg = 'Check cannot be performed since frequency step is not constant'
   MSG_WARNING(msg)
   RETURN
  end if 
 end do

!Check whether omega(1) is small or not
 if (omega(1) > 0.1/Ha_eV) then
  msg='Check cannot be performed since first frequency on the grid > 0.1 eV'
  MSG_WARNING(msg)
  RETURN 
 end if 

!If eps(nomega) is not 0 warn
 if (eps(nomega) > 0.1) then
  write(msg,'(a,f8.4,3a,f8.2,2a)')&
&  '  Im epsilon for omega = ',omega(nomega)*Ha_eV,' eV',ch10,&
&  '  is not yet zero, epsilon_2 = ',eps(nomega),ch10,&
&  '  f-sum rule test could give wrong results'
  MSG_WARNING(msg)
 end if

!Check f-sum rule using naive integration
 SELECT CASE (method)

 CASE (0)
  acc=zero
  do ip=1,nomega
   wwp=omega(ip)
   acc=acc+wwp*e2(ip)
  end do
  acc=domega*acc
  
  omegaplasmaeff=SQRT(acc*2/pi)

! Perform Kramers-Kronig using Simpson integration    
! Simpson O(1/N^4), from NumRec in C p 134  NumRec in Fortran p 128
 CASE (1)
  do ip=1,nomega
   wwp=omega(ip)
   fs(ip)=wwp*e2(ip) 
  end do

  call simspon_int(nomega,domega,fsint)
  omegaplasmaeff=SQRT(fsint(nomega)*2/pi)

 CASE DEFAULT 
  write(msg,'(a,i3)')' wrong value for method ',method
  MSG_BUG(msg)
 END SELECT

 check=ABS((omegaplasmaeff-omegaplasma))/omegaplasma

!Write data
 write(msg,'(2a,f6.2,3a,f6.2,3a,f6.2,2a)'),ch10,&
& '( omega_plasma     = ',omegaplasma*Ha_eV,  ' [eV] )',ch10,&
& '( omega_plasma^eff = 'omegaplasmaeff*Ha_eV,' [eV] )',ch10,&
& '( the f-sum rule is verified within ',check*100,'% )',ch10
 call wrtout(std_out,msg,'COLL')

end subroutine fsumrule
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/check_kramerskronig
!! NAME
!!  check_kramerskronig
!!
!! FUNCTION
!!   check Kramers Kronig
!!   \int_0^\infty d\omega' frac{\omega'}{\omega'^2 - \omega^2} 
!!   Im \epsilon(\omega') = Re \epsilon(\omega)
!!
!! COPYRIGHT
!! Copyright (C) 2009-2016 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  n=Number of frequency points.
!!  eps(n)=Dielectric function.
!!  o(n)=Frequency mesh.
!!
!! OUTPUT
!!  Only checking.
!!
!! PARENTS
!!      m_exc_spectra
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine check_kramerskronig(n,o,eps)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_numeric_tools,  only : simpson_cplx

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'check_kramerskronig'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
!arrays
 complex(dpc),intent(in) :: eps(n)
 real(dp),intent(in) :: o(n)

!Local variables ------------------------------
!scalars
 integer :: ii,ip
 real(dp) :: omega,omegap,domega,kk,kkrms,eav
 complex(dpc) c
 character(len=500) :: msg
!arrays
 real(dp) :: e1kk(n)
 complex(dpc) :: intg(n)

!************************************************************************
! init jmb
 e1kk=zero
 intg=(zero,zero)

! calculate domega step and verify all  
 domega = (o(n) - o(1)) / (n-1)

 do ii=2,n
  if (domega-(o(ii)-o(ii-1)) > tol3) then 
    MSG_WARNING("Frequency mesh not linear. Returning")
    return
  end if
 end do

 if(o(1) > 0.1/Ha_eV) then
   MSG_WARNING("First frequency is not zero. Returning")
   return
 end if

 if (aimag(eps(n)) > 0.1) then
   write(msg,'(a,f12.6,3a,f12.6,2a)')&
&   ' Im epsilon for omega= ',o(n)*Ha_eV,'eV',ch10,&
&   ' is not yet zero, epsilon_2= ',aimag(eps(n)),ch10,&
&   ' Kramers Kronig test could give wrong results. '
   MSG_WARNING(msg)
 end if
      
! Fill array for kramers kronig.     
 do ii=1,n
   omega=o(ii)
   c = (0.0,0.0)
   do ip=1,n
     if(ip == ii) cycle
     omegap = o(ip)
     c = c + omegap / (omegap**2-omega**2) * aimag(eps(ip))
   end do
   e1kk(ii) = one + two/pi * domega*real(c)
 end do
      
!perform kramers kronig with simpson integration    
 do ii=1,n
   omega=o(ii)
   do ip=1,n
     if (ip==ii) cycle
     omegap = o(ip)
     intg(ip) = omegap / (omegap**2 - omega**2) * aimag(eps(ip))
   end do
   c = simpson_cplx(n,domega,intg)
   e1kk(ii) = one + two/pi * real(c)
 end do

!verify kramers kronig
 eav=zero; kk=zero; kkrms=zero
 do ii=1,n
   kk = kk + abs(real(eps(ii)) - e1kk(ii))
   kkrms = kkrms +(real(eps(ii)) - e1kk(ii))*(real(eps(ii)) - e1kk(ii))
   eav = eav + abs(real(eps(ii)))
 end do

 eav = eav/n
 kk = (kk/n)/eav
 kkrms = (kkrms/n) / (eav*eav)

 kk = abs(real(eps(1)) - e1kk(1)) / real(eps(1))

! write data
 write(msg,'(a,f7.2,a)')" The Kramers-Kronig is verified within ",100*kk,"%"
 call wrtout(std_out,msg,"COLL")

! write(std_out,'("# Kramers Kronig calculation of epsilon1")')
! write(std_out,'("# omega   epsilon1  epsilon1kk")')
! do ii=1,n
!   write(std_out,'(f7.3,2e15.7)') o(ii)*Ha_eV, real(eps(ii)), e1kk(ii)
! end do
      
end subroutine check_kramerskronig
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/check_fsumrule
!! NAME
!!  check_fsumrule
!!
!! FUNCTION
!!   check f-sum rule
!!   \int_0^\infty d\omega \omega Im \epsilon_GG'(q,\omega) =
!!   = \frac{1}{2} \pi \omega_p^2  \frac{\rho(G-G')}{\rho(0)} 
!!   versor(q+G) \dot versor(q+G')
!!   for q = G = G' = 0, it reads:
!!   \int_0^\infty d\omega \omega Im \epsilon_00(q=0,\omega) = 
!!   = \pi \omega_p^2 / 2
!!   calculate only the second one
!!   calculate the integral to evaluate an omega_plasma^eff to compare with omega_plasma
!!
!! COPYRIGHT
!! Copyright (C) 2009-2016 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  n=Number of frequencies.
!!  o(n)=Frequency mesh.
!!  e2(n)=imaginary part of epsilon_00
!!  omegaplasma=Drude plasma frequency.
!!
!! OUTPUT
!!  Only checking.
!!
!! PARENTS
!!      m_exc_spectra
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine check_fsumrule(n,o,e2,omegaplasma)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_numeric_tools,  only : simpson_cplx

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'check_fsumrule'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n 
 real(dp),intent(in) :: omegaplasma
!arrays
 real(dp),intent(in) :: o(n),e2(n) 

!Local variables ------------------------------
!scalars
 integer :: ii,ip
 real(dp) :: omegap,domega,integral,omegaplasmaeff,fsumrule
 character(len=500) :: msg
!arrays
 complex(dpc) :: intg(n)

!************************************************************************

! calculate domega step and verify      
 domega = (o(n) - o(1)) / (n-1)

 do ii=2,n
   if (domega-(o(ii)-o(ii-1)) > tol3) then
     MSG_WARNING("Frequency mesh not linear. Returning")
     return
   end if
 end do

 if (o(1) > 0.1/Ha_eV) then
   MSG_WARNING("First frequency is not zero. Returning")
   return                                                
 end if
 
 if (e2(n) > 0.1) then
   write(msg,'(a,f12.6,3a,f12.6,2a)')&
&   ' Im epsilon for omega= ',o(n)*Ha_eV,' eV ',ch10,&
&   ' is not yet zero, epsilon_2= ',e2(n),ch10,&
&   ' f-sum rule test could give wrong results.'
   MSG_WARNING(msg)
 end if
      
! integrate to obtain f-sum rule
 integral=zero
 do ip=1,n
   omegap=o(ip)
   integral = integral + omegap * e2(ip)
 end do
 integral = domega * integral
      
!integrate with simpson to obtain f-sum rule   
 do ip = 1, n
   omegap = o(ip)
   intg(ip) = omegap * e2(ip) 
 end do

 integral = real(simpson_cplx(n,domega,intg))
 if(integral < 0) then
   MSG_ERROR("The integral of the imaginary of dielectric function is negative !!!")
 else
   omegaplasmaeff = sqrt(integral*two/pi)
 end if

 fsumrule = abs((omegaplasmaeff - omegaplasma)) / omegaplasma

! write data
 write(msg,'(3(a,f6.2,2a))')&
&  " omega_plasma     = ",omegaplasma*Ha_eV,   " [eV]",ch10,&
&  " omega_plasma^eff = ",omegaplasmaeff*Ha_eV," [eV]",ch10,&
&  " the f-sum rule is verified within ",fsumrule*100,"%",ch10
 call wrtout(std_out,msg,"COLL")

end subroutine check_fsumrule
!!***


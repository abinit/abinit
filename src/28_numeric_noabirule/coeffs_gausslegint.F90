!{\src2tex{textfont=tt}}
!!****f* ABINIT/coeffs_gausslegint
!! NAME
!!  coeffs_gausslegint
!!
!! FUNCTION
!! Compute the coefficients (supports and weights)
!! for Gauss-Legendre integration.
!! Inspired by a routine due to G. Rybicki.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2016 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! xmin=lower bound of integration
!! xmax=upper bound of integration
!! n=order of integration
!!
!! OUTPUT
!! x(n)=array of support points
!! weights(n)=array of integration weights
!!
!! PARENTS
!!      calc_rpa_functional,calc_sigc_me,integrho,integvol,m_numeric_tools
!!      screening,surf
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine coeffs_gausslegint(xmin,xmax,x,weights,n)

 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'coeffs_gausslegint'
!End of the abilint section

      implicit none
      
      integer :: n,i,j 
      double precision :: tol,xmin,xmax,x(n),weights(n)
      double precision :: z,z1,pi,xmean,p1,p2,p3,pp,xl

      tol=1.d-13
      pi=4.d0*atan(1.d0)

      xl=(xmax-xmin)*0.5d0
      xmean=(xmax+xmin)*0.5d0

      do i=1,(n+1)/2
       z=cos(pi*(i-0.25d0)/(n+0.5d0))
      
       do 

         p1=1.d0
         p2=0.d0
      
         do j=1,n
          
          p3=p2
          p2=p1
          p1=((2.d0*j - 1.d0)*z*p2 - (j-1.d0)*p3)/j
        
         enddo
       
         pp=n*(p2-z*p1)/(1.0d0-z**2)
         z1=z
         z=z1-p1/pp
         
         if(abs(z-z1) < tol) exit

       enddo

       x(i)=xmean-xl*z
       x(n+1-i)=xmean+xl*z
       weights(i)=2.d0*xl/((1.d0-z**2)*pp**2)
       weights(n+1-i)=weights(i)

      enddo

end subroutine
!!***

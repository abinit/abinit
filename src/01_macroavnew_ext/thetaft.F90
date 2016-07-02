#if defined HAVE_CONFIG_H
#include "config.h"
#endif

        subroutine thetaft(n,L,lav,ft)
 
! Copyright (C) 1999-2003 (P. Ordejon, J. Junquera)
! This file is distributed under the terms of the
! GNU General Public License, see ~abinit/COPYING
! or http://www.gnu.org/copyleft/gpl.txt .
 
!
! Calculates the Fourier coefficients of the convoluting function
!    w(x) = (1/lav) theta(lav/2 - abs(x))
! asuming it is periodic with period L:
!
!       |                                          |
!  1/lav|__________                      __________|
!       |          |                     |         |
!       |__________|_____________________|_________|___
!       0         lav/2               L-lav/2     L
!

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'thetaft'
!End of the abilint section

        implicit none

        integer n,j
        real*8 L,lav
        real*8 ft(2*n)
        real*8 pi

        pi=4.0d0*datan(1.0d0)

        ft(1)=n/L
        ft(2)=0.0

        do j=2,n/2+1
          ft(2*(j-1)+1) = (n/(pi*lav*(j-1)/2.))*&
     &                    sin(pi*lav*(j-1)/2./L)*&
     &                    cos(pi*(L-lav/2.)*(j-1)/L)*&
     &                    cos(pi*L*(j-1)/L)
          ft(2*(j-1)+2) = (n/(pi*lav*(j-1)/2.))*&
     &                    sin(pi*lav*(j-1)/2./L)*&
     &                    cos(pi*(L-lav/2.)*(j-1)/L)*&
     &                    sin(pi*L*(j-1)/L)
        enddo

        do j=n/2+2,n
          ft(2*(j-1)+1) = (n/(pi*lav*(-n+j-1)/2.))*&
     &                    sin(pi*lav*(-n+j-1)/2./L)*&
     &                    cos(pi*(L-lav/2.)*(-n+j-1)/L)*&
     &                    cos(pi*L*(-n+j-1)/L)
          ft(2*(j-1)+2) = (n/(pi*lav*(-n+j-1)/2.))*&
     &                    sin(pi*lav*(-n+j-1)/2./L)*&
     &                    cos(pi*(L-lav/2.)*(-n+j-1)/L)*&
     &                    sin(pi*L*(-n+j-1)/L)
        enddo


        return
        end

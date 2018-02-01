!{\src2tex{textfont=tt}}
!!****f* ABINIT/jacobi
!! NAME
!!  jacobi
!!
!! FUNCTION
!!     Computes all eigenvalues and eigenvectors of a real symmetric matrix a,
!!     which is of size n by n, stored in a physical np by np array. On output,
!!     elements of a above the diagonal are destroyed. d returns the
!!     eigenvalues of a in its first n elements. v is a matrix with the same
!!     logical and physical dimensions as a, whose columns contain, on output,
!!     the normalized eigenvectors of a. nrot returns the number of Jacobi
!!     rotations that were required.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2018 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      conducti_nc,critic
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine jacobi(a,n,np,d,v,nrot)

 use defs_basis, only : std_out

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'jacobi'
!End of the abilint section

 implicit none
!Arguments
 integer :: n,np,nrot
 real*8 :: a(np,np),d(np),v(np,np)
!Local variables
 integer, parameter :: NMAX=500
 integer i,ip,iq,j
 real*8 c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
 do ip=1,n 
   do iq=1,n
     v(ip,iq)=0.
   enddo
   v(ip,ip)=1.
 enddo
 do ip=1,n
   b(ip)=a(ip,ip) 
   d(ip)=b(ip)
   z(ip)=0. 
 enddo
 nrot=0
 do i=1,50
   sm=0.
   do ip=1,n-1 
     do iq=ip+1,n
       sm=sm+abs(a(ip,iq))
     enddo
   enddo
   if(sm.eq.0.)return 
   if(i.lt.4)then
     tresh=0.2*sm/n**2 
   else
     tresh=0. 
   endif
   do ip=1,n-1
     do iq=ip+1,n
       g=100.*abs(a(ip,iq))
       if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))) &
&          .and.(abs(d(iq))+g.eq.abs(d(iq))))then
            a(ip,iq)=0.
       else if(abs(a(ip,iq)).gt.tresh)then
         h=d(iq)-d(ip)
         if(abs(h)+g.eq.abs(h))then
           t=a(ip,iq)/h 
         else
           theta=0.5*h/a(ip,iq) 
           t=1./(abs(theta)+sqrt(1.+theta**2))
           if(theta.lt.0.)t=-t
         endif
         c=1./sqrt(1+t**2)
         s=t*c
         tau=s/(1.+c)
         h=t*a(ip,iq)
         z(ip)=z(ip)-h
         z(iq)=z(iq)+h
         d(ip)=d(ip)-h
         d(iq)=d(iq)+h
         a(ip,iq)=0.
         do j=1,ip-1
           g=a(j,ip)
           h=a(j,iq)
           a(j,ip)=g-s*(h+g*tau)
           a(j,iq)=h+s*(g-h*tau)
         enddo
         do j=ip+1,iq-1
           g=a(ip,j)
           h=a(j,iq)
           a(ip,j)=g-s*(h+g*tau)
           a(j,iq)=h+s*(g-h*tau)
         enddo
         do j=iq+1,n 
           g=a(ip,j)
           h=a(iq,j)
           a(ip,j)=g-s*(h+g*tau)
           a(iq,j)=h+s*(g-h*tau)
         enddo
         do j=1,n
           g=v(j,ip)
           h=v(j,iq)
           v(j,ip)=g-s*(h+g*tau)
           v(j,iq)=h+s*(g-h*tau)
         enddo
         nrot=nrot+1
       endif
     enddo
   enddo
   do ip=1,n
     b(ip)=b(ip)+z(ip)
     d(ip)=b(ip) 
     z(ip)=0. 
   enddo
 enddo
 write(std_out,*) 'too many iterations in jacobi'

end subroutine jacobi
!!***

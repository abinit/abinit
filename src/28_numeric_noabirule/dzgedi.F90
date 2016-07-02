!{\src2tex{textfont=tt}}
!!****f* ABINIT/dzegdi
!! NAME
!!  dzgedi
!!
!! FUNCTION
!!  This routine is the clone of zgefa.F90 using real*8 a(2) instead of complex*16
!!  for the purpose of ABINIT
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2016 ABINIT group (XG)
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
!!      berryphase,dfptnl_mv,qmatrix,relaxpol,smatrix,uderiv
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dzgedi(a,lda,n,ipvt,det,work,job)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dzgedi'
!End of the abilint section

      implicit none

      integer :: lda,n,ipvt(n),job
      real*8 :: a(2,lda,n),det(2,2),work(2,n)
!
!     zgedi computes the determinant and inverse of a matrix
!     using the factors computed by zgeco or zgefa.
!
!     on entry
!
!        a       complex*16(lda, n)
!                the output from zgeco or zgefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        ipvt    integer(n)
!                the pivot vector from zgeco or zgefa.
!
!        work    complex*16(n)
!                work vector.  contents destroyed.
!
!        job     integer
!                = 11   both determinant and inverse.
!                = 01   inverse only.
!                = 10   determinant only.
!
!     on return
!
!        a       inverse of original matrix if requested.
!                otherwise unchanged.
!
!        det     complex*16(2)
!                determinant of original matrix if requested.
!                otherwise not referenced.
!                determinant = det(1) * 10.0**det(2)
!                with  1.0 .le. cabs1(det(1)) .lt. 10.0
!                or  det(1) .eq. 0.0 .
!
!     error condition
!
!        a division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        it will not occur if the subroutines are called correctly
!        and if zgeco has set rcond .gt. 0.0 or zgefa has set
!        info .eq. 0 .
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     internal variables
!
      double precision :: r(2),rk(2),rkj(2)
      double precision :: ten,rinv2,rabs
      integer :: i,j,k,kb,kp1,l,nm1
!
!     compute determinant
!
      if (job/10 .eq. 0) go to 70
         det(1,1) = 1.0d0; det(2,1) = 0.0d0
         det(1,2) = 0.0d0; det(2,2) = 0.0d0
         ten = 10.0d0
         do i = 1, n
            if (ipvt(i) .ne. i) then
                det(1,1) = -det(1,1)
                det(2,1) = -det(2,1)
            end if
            r(1)=det(1,1); r(2)=det(2,1)
            det(1,1) = r(1)*a(1,i,i)-r(2)*a(2,i,i)
            det(2,1) = r(2)*a(1,i,i)+r(1)*a(2,i,i)
!        ...exit
            rabs = abs(det(1,1))+abs(det(2,1))
            if (rabs .eq. 0.0d0) go to 60
   10       continue
            rabs = abs(det(1,1))+abs(det(2,1))
            if (rabs .ge. 1.0d0) go to 20
               det(1,1) = ten*det(1,1); det(2,1) = ten*det(2,1)
               det(1,2) = det(1,2) - 1.0d0
            go to 10
   20       continue
   30       continue
            rabs = abs(det(1,1))+abs(det(2,1))
            if (rabs .lt. ten) go to 40
               det(1,1) = det(1,1)/ten; det(2,1) = det(2,1)/ten
               det(1,2) = det(1,2) + 1.0d0
            go to 30
   40       continue
         end do
   60    continue
   70 continue
!
!     compute inverse(u)
!
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            !a(k,k) = (1.0d0,0.0d0)/a(k,k)
            !t = -a(k,k)
            !call zscal(k-1,t,a(1,k),1)
            rinv2 = 1.d0/(a(1,k,k)**2+a(2,k,k)**2)
            a(1,k,k) =  rinv2*a(1,k,k)
            a(2,k,k) = -rinv2*a(2,k,k)
            rk(1) = -a(1,k,k); rk(2) = -a(2,k,k)
            do i=1,k-1
               r(1)=a(1,i,k)
               r(2)=a(2,i,k)
               a(1,i,k)=rk(1)*r(1)-rk(2)*r(2)
               a(2,i,k)=rk(1)*r(2)+rk(2)*r(1)
            end do
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               !t = a(k,j)
               !a(k,j) = (0.0d0,0.0d0)
               !call zaxpy(k,t,a(1,k),1,a(1,j),1)
               rkj(1) = a(1,k,j); rkj(2) = a(2,k,j)
               a(1,k,j) = 0.d0; a(2,k,j) = 0.d0
               do i=1,k
                  a(1,i,j)=rkj(1)*a(1,i,k)-rkj(2)*a(2,i,k)+a(1,i,j)
                  a(2,i,j)=rkj(2)*a(1,i,k)+rkj(1)*a(2,i,k)+a(2,i,j)
               end do
   80       continue
   90       continue
  100    continue
  do i=1,n
  end do
!
!        form inverse(u)*inverse(l)
!
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(1,i) = a(1,i,k); work(2,i) = a(2,i,k)
               a(1,i,k) = 0.0d0; a(2,i,k) = 0.d0
  110       continue
            do 120 j = kp1, n
               r(1) = work(1,j); r(2) = work(2,j)
               !call zaxpy(n,t,a(1,j),1,a(1,k),1)
               do i=1,n
                  a(1,i,k)=r(1)*a(1,i,j)-r(2)*a(2,i,j)+a(1,i,k)
                  a(2,i,k)=r(2)*a(1,i,j)+r(1)*a(2,i,j)+a(2,i,k)
               end do
  120       continue
            l = ipvt(k)
            if (l .ne. k) then
               !call zswap(n,a(1,k),1,a(1,l),1)
               do i=1,n
                  r(1) = a(1,i,k); r(2) = a(2,i,k)
                  a(1,i,k) = a(1,i,l); a(2,i,k) = a(2,i,l)
                  a(1,i,l) = r(1); a(2,i,l) = r(2)
               end do
            end if
  130    continue
  140    continue
  150 continue

end subroutine dzgedi
!!***

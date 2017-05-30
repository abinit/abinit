!{\src2tex{textfont=tt}}
!!****f* ABINIT/dzgefa
!! NAME
!!  dzgefa
!!
!! FUNCTION
!!   This routine is the clone of zgefa.F90 using real*8 a(2) instead of complex*16
!!   for the purpose of ABINIT (2008,TD)
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2017 ABINIT group (XG)
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

subroutine dzgefa(a,lda,n,ipvt,info)

 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dzgefa'
!End of the abilint section

 implicit none

!Arguments
 integer :: lda,n,ipvt(n),info
 real*8  :: a(2,lda,n)
!
!     zgefa factors a complex*16 matrix by gaussian elimination.
!
!     dzgefa is usually called by zgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for zgeco) = (1 + 9/n)*(time for zgefa) .
!
!     on entry
!
!        a       complex*16(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that zgesl or zgedi will divide by zero
!                     if called.  use  rcond  in zgeco for a reliable
!                     indication of singularity.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     internal variables
!
!Local variables
 real*8 :: r(2),rk(2),rlj(2)
 real*8 :: rinv2,rmax,rabs
 integer :: i,j,k,kp1,l,nm1

!
!     gaussian elimination with partial pivoting
!
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
!
!        find l = pivot index
!
         !l = izamax(n-k+1,a(k,k),1) + k - 1
         rmax=0.d0
         l=0
         do i=k,n
            rabs=abs(a(1,i,k))+abs(a(2,i,k))
            if(rmax<=rabs) then
              rmax=rabs
              l=i
            end if
         end do
         ipvt(k) = l
!
!        zero pivot implies this column already triangularized
!
         if (abs(a(1,l,k))+abs(a(2,l,k)) .eq. 0.0d0) go to 40
!
!           interchange if necessary
!
            if (l .eq. k) go to 10
               r(1) = a(1,l,k); r(2) = a(2,l,k)
               a(1,l,k) = a(1,k,k); a(2,l,k) = a(2,k,k)
               a(1,k,k) = r(1); a(2,k,k) = r(2)
   10       continue
!
!           compute multipliers
!
            rinv2 = 1.d0/(a(1,k,k)**2+a(2,k,k)**2)
            rk(1) = -rinv2*a(1,k,k)
            rk(2) =  rinv2*a(2,k,k)
            !call zscal(n-k,t,a(k+1,k),1)
            do i=k+1,n
               r(1)=a(1,i,k)
               r(2)=a(2,i,k)
               a(1,i,k)=rk(1)*r(1)-rk(2)*r(2)
               a(2,i,k)=rk(1)*r(2)+rk(2)*r(1)
            end do
!
!           row elimination with column indexing
!
            do j = kp1, n
               rlj(1) = a(1,l,j); rlj(2) = a(2,l,j)
               if (l .eq. k) go to 20
                  a(1,l,j) = a(1,k,j); a(2,l,j) = a(2,k,j)
                  a(1,k,j) = rlj(1); a(2,k,j) = rlj(2)
   20          continue
               !call zaxpy(n-k,t,a(1,k+1,k),1,a(1,k+1,j),1)
               do i=k+1,n
                  a(1,i,j)=rlj(1)*a(1,i,k)-rlj(2)*a(2,i,k)+a(1,i,j)
                  a(2,i,j)=rlj(2)*a(1,i,k)+rlj(1)*a(2,i,k)+a(2,i,j)
               end do
            end do
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (abs(a(1,n,n))+abs(a(2,n,n)) .eq. 0.0d0) info = n

end subroutine dzgefa
!!***

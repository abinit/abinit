!{\src2tex{textfont=tt}}
!!****f* ABINIT/ludcmp
!! NAME
!!  ludcmp
!!
!! FUNCTION
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
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


SUBROUTINE ludcmp(a,n,np,indx,id,info)

 use defs_basis, only : std_out

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ludcmp'
!End of the abilint section

      implicit none

      INTEGER n,np,indx(n),NMAX,id,info
      REAL*8 a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)

!     Given a matrix a(1:n,1:n), with physical dimension np by np, this
!     routine replaces it by the LU decomposition of a rowwise permutation of
!     itself. a and n are input. a is output, arranged as in equation (2.3.14)
!     above; indx(1:n) is an output vector that records the row permutation
!     effected by the partial pivoting; id is output as +- 1 depending on
!     whether the number of row interchanges was even or odd,
!     respectively. This routine is used in combination with lubksb to solve
!     linear equations or invert a matrix.

      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX) 

!      write(std_out,*) 'ENTERING LUDCMP...'
!      write(std_out,*) 'in ludcmp n=',n,' np=',np
!      write(std_out,201) ((a(i,j),j=1,n),i=1,n)
! 201  FORMAT('A in ludcmp ',/,3F16.8,/,3F16.8,/,3F16.8)
      id=1
      info=0
      do i=1,n 
        aamax=0.
        do j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo
        if (aamax.eq.0.) then
          write(std_out,*) 'LUDCMP: singular matrix !!!'
          do j=1,3
            write(std_out,*) (a(j,k),k=1,3)
          enddo
          info=1
          return
!          stop 'singular matrix in ludcmp' 
        endif
        vv(i)=1./aamax 
      enddo 
      do j=1,n 
        do i=1,j-1 
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo 
          a(i,j)=sum
        enddo 
        aamax=0. 
        do i=j,n 
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo 
          a(i,j)=sum
          dum=vv(i)*abs(sum) 
          if (dum.ge.aamax) then 
            imax=i
            aamax=dum
          endif
        enddo 
        if (j.ne.imax)then 
          do  k=1,n 
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo 
          id=-id 
          vv(imax)=vv(j) 
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then 
          dum=1./a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo 
!      write(std_out,*) 'LEAVING LUDCMP...'
      return
END
!!***

!!****f* ABINIT/lubksb
!! NAME
!!  lubksb
!!
!! FUNCTION
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
!!
!! CHILDREN
!!
!! SOURCE


SUBROUTINE lubksb(a,n,np,indx,b)

 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lubksb'
!End of the abilint section

      IMPLICIT NONE

      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)

!     Solves the set of n linear equations A . X = B. Here a is input, not as
!     the matrix A but rather as its LU decomposition, determined by the
!     routine ludcmp. indx is input as the permutation vector returned by
!     ludcmp. b(1:n) is input as the right-hand side vector B, and returns
!     with the solution vector X. a, n, np, and indx are not modified by this
!     routine and can be left in place for successive calls with different
!     right-hand sides b. This routine takes into account the possibility that
!     b will begin with many zero elements, so it is efficient for use in
!     matrix inversion.

      INTEGER i,ii,j,ll
      REAL*8 sum
!      write(std_out,*) 'ENTERING LUBKSB...'
!      write(std_out,201) ((a(i,j),j=1,n),i=1,n)
! 201  FORMAT('A in lubksb ',/,3F16.8,/,3F16.8,/,3F16.8)

      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo 
        else if (sum.ne.0.) then
          ii=i 
        endif
        b(i)=sum
      enddo 
      do i=n,1,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i) 
      enddo
!      write(std_out,*) 'LEAVING LUBKSB...'
      return 
END SUBROUTINE LUBKSB
!!***


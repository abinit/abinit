!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_kpt_fullbz
!! NAME
!! get_kpt_fullbz
!!
!! FUNCTION
!! create full grid of kpoints from kptrlatt and shiftk
!!
!! COPYRIGHT
!! Copyright (C) 2002-2017 ABINIT group (MVer,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  kptrlatt(3,3)=lattice vectors for full kpoint grid
!!  nkpt_fullbz=number of kpoints in full brillouin zone
!!  nshiftk=number of kpoint grid shifts
!!  shiftk(3,nshiftk)=kpoint shifts
!!
!! OUTPUT
!!  kpt_fullbz(3,nkpt_fullbz)=kpoints in full brillouin zone
!!
!! NOTES
!!
!! PARENTS
!!      get_full_kgrid,invars2
!!
!! CHILDREN
!!      mati3det,matr3inv,wrap2_pmhalf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine get_kpt_fullbz(kpt_fullbz,kptrlatt,nkpt_fullbz,nshiftk,shiftk)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_numeric_tools,   only : wrap2_pmhalf

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_kpt_fullbz'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt_fullbz,nshiftk
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk)
 real(dp),intent(out) :: kpt_fullbz(3,nkpt_fullbz)

!Local variables-------------------------------
!scalars
 integer, parameter :: max_number_of_prime=47,mshiftk=210
 integer :: det,ii,ikshft,iprim,jj,kk,nn
 character(len=500) :: message

!arrays
 integer :: boundmax(3),boundmin(3),common_factor(3)
 integer, parameter :: prime_factor(max_number_of_prime)=(/2,3,5,7,9, 11,13,17,19,23,&
&  29,31,37,41,43, 47,53,59,61,67,&
&  71,73,79,83,89, 97,101,103,107,109,&
&  113,127,131,137,139, 149,151,157,163,167,&
&  173,179,181,191,193, 197,199/)
 real(dp) :: k1(3),k2(3),klatt(3,3),rlatt(3,3),shift(3),test_rlatt(3,3)

! *********************************************************************

!Identify first factors that can be used to rescale the three kptrlatt vectors
!Only test a large set of prime factors, though ...
 do jj=1,3
   common_factor(jj)=1
   rlatt(:,jj)=kptrlatt(:,jj)
   do iprim=1,max_number_of_prime
     test_rlatt(:,jj)=rlatt(:,jj)/dble(prime_factor(iprim))
!    If one of the components is lower than 1 in absolute value, then it is not worth to continue the search.
     if(minval(abs(abs(test_rlatt(:,jj))-half))<half-tol8)exit
     do
       if(sum(abs(test_rlatt(:,jj)-nint(test_rlatt(:,jj)) ))<tol8)then
         common_factor(jj)=prime_factor(iprim)*common_factor(jj)
         rlatt(:,jj)=rlatt(:,jj)/dble(prime_factor(iprim))
         test_rlatt(:,jj)=test_rlatt(:,jj)/dble(prime_factor(iprim))
       else
         exit
       end if  
     end do
   end do
 end do
 call mati3det(kptrlatt,det)
 det=det/(common_factor(1)*common_factor(2)*common_factor(3))

 rlatt(:,:)=kptrlatt(:,:)
 call matr3inv(rlatt,klatt)
!Now, klatt contains the three primitive vectors of the k lattice,
!in reduced coordinates. One builds all k vectors that
!are contained in the first Brillouin zone, with coordinates
!in the interval [0,1[ . First generate boundaries of a big box.
!In order to generate all possible vectors in the reciprocal space, 
!one must consider all multiples of the primitive ones, until a vector with only integers is found.
!The maximum bound is the scale of the corresponding kptrlatt vector, times the determinant of kptrlatt. Also consider negative vectors.
!On this basis, compute the bounds.
 do jj=1,3
!  To accomodate the shifts, boundmin starts from -1
!  Well, this is not a complete solution ...
   boundmin(jj)=-1-common_factor(jj)*abs(det)
   boundmax(jj)=common_factor(jj)*abs(det)
 end do

 nn=1
 do kk=boundmin(3),boundmax(3)
   do jj=boundmin(2),boundmax(2)
     do ii=boundmin(1),boundmax(1)
       do ikshft=1,nshiftk

!        Coordinates of the trial k point with respect to the k primitive lattice
         k1(1)=ii+shiftk(1,ikshft)
         k1(2)=jj+shiftk(2,ikshft)
         k1(3)=kk+shiftk(3,ikshft)

!        Reduced coordinates of the trial k point
         k2(:)=k1(1)*klatt(:,1)+k1(2)*klatt(:,2)+k1(3)*klatt(:,3)

!        Eliminate the point if outside [0,1[
         if(k2(1)<-tol10)cycle ; if(k2(1)>one-tol10)cycle
         if(k2(2)<-tol10)cycle ; if(k2(2)>one-tol10)cycle
         if(k2(3)<-tol10)cycle ; if(k2(3)>one-tol10)cycle

!        Wrap the trial values in the interval ]-1/2,1/2] .
         call wrap2_pmhalf(k2(1),k1(1),shift(1))
         call wrap2_pmhalf(k2(2),k1(2),shift(2))
         call wrap2_pmhalf(k2(3),k1(3),shift(3))
         if(nn > nkpt_fullbz) then
           write (message,'(a,i0)')' nkpt_fullbz mis-estimated, exceed nn=',nn
           MSG_BUG(message)
         end if
         kpt_fullbz(:,nn)=k1(:)
         nn=nn+1
       end do
     end do
   end do
 end do
 nn = nn-1

 if (nn /= nkpt_fullbz) then
   write (message,'(2(a,i0))')' nkpt_fullbz= ',nkpt_fullbz,' underestimated  nn=',nn
   MSG_BUG(message)
 end if

end subroutine get_kpt_fullbz
!!***

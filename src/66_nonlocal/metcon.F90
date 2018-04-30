!{\src2tex{textfont=tt}}
!!****f* ABINIT/metcon
!! NAME
!! metcon
!!
!! FUNCTION
!! Carries out specialized metric tensor contractions needed for
!! l=0,1,2,3 nonlocal Kleinman-Bylander pseudopotential operation.
!! Full advantage is taken of the full permutational symmetry of these
!! tensors.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  rank=0,1,2, or 3 = rank of input tensor aa
!!  gmet(3,3)=metric tensor (array is symmetric but stored as 3x3)
!!  aa(2,(rank+1)*(rank+2)/2)=unique elements of complex input tensor
!!
!! OUTPUT
!!  bb(2,(rank+1)*(rank+2)/2)=unique elements of complex output tensor
!!
!! NOTES
!! All tensors are stored in a compressed storage mode defined below;
!! input and output conform to this scheme.
!! When tensor elements occur repeatedly due to symmetry, the
!! WEIGHT IS INCLUDED in the output tensor element to simplify later
!! contractions with other tensors of the same rank and form, i.e. the
!! next contraction is then simply a dot product over the unique elements.
!!
!! Definitions of the contractions:
!!
!! rank=0:  bb = aa  (simply copy the scalar value, no use of gmet)
!!
!! rank=1:  bb(i)= $gmet(i,l) aa(l)$  (3 elements in, 3 elements out)
!!
!! rank=2:  bb(i,j)= $[{3 \over 2} gmet(i,l) gmet(j,m) - {1 \over 2} gmet(i,j) gmet(l,m)] aa(l,m)$
!!          (6 elements in, 6 elements out)
!!
!! rank=3:  bb(i,j,k)= $[{5 \over 2} g(i,l) g(j,m) g(k,n) - {3 \over 2} g(i,j) g(l,m) g(k,n)] aa(l,m,n)$
!!          (10 elements in, 10 elements out)
!!         In this rank 3 case, the second term is NOT symmetric in all
!!         permutations of i,j,k, but the final tensor b(ijk) may be
!!         symmetrized over all permutations because it will be
!!         contracted with a completely symmetric tensor.
!!
!! The compressed storage scheme is based on storing a symmetric 3x3 matrix as
!!     (1 . .)
!!     (6 2 .)
!!     (5 4 3)
!! which leads to the following mappings for all ranks
!! where the compressed storage index is to the right of the arrow:
!! rank=0 1->1 (only a scalar)
!! rank=1 1->1 2->2 3->3 (standard vector, no compression)
!! rank=2 11->1 22->2 33->3 32->4 31->5 21->6
!!  weights   1     1     1     2     2     2
!! rank=3 111->1 221->2 331->3 321->4 311->5 211->6 222->7 332->8 322->9 333->10
!!  weights    1      3      3      6      3      3      1      3      3       1
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine metcon(rank,gmet,aa,bb)

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'metcon'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rank
!arrays
 real(dp),intent(in) :: aa(2,((rank+1)*(rank+2))/2),gmet(3,3)
 real(dp),intent(out) :: bb(2,((rank+1)*(rank+2))/2)

!Local variables-------------------------------
!scalars
 integer :: ii,jj
 real(dp) :: scalar,tmpiii,tmpijk
 character(len=500) :: message
!arrays
 real(dp) :: vector(3)

! *************************************************************************

!This statement function defines the l=3 contraction in
!terms of the free indices of the contracted tensor (re and im)
! coniii(ii,i1,i2,i3)=gmet(i1,1)*gmet(i2,1)*gmet(i3,1)*aa(ii,1)+&
!& gmet(i1,2)*gmet(i2,2)*gmet(i3,2)*aa(ii,7)+&
!& gmet(i1,3)*gmet(i2,3)*gmet(i3,3)*aa(ii,10)
! conijk(ii,i1,i2,i3)=aa(ii,4)*&
!& (gmet(i1,1)*gmet(i2,2)*gmet(i3,3)+&
!& gmet(i1,2)*gmet(i2,3)*gmet(i3,1)+&
!& gmet(i1,3)*gmet(i2,1)*gmet(i3,2)+&
!& gmet(i1,3)*gmet(i2,2)*gmet(i3,1)+&
!& gmet(i1,1)*gmet(i2,3)*gmet(i3,2)+&
!& gmet(i1,2)*gmet(i2,1)*gmet(i3,3))
! con(ii,i1,i2,i3)=coniii(ii,i1,i2,i3)+conijk(ii,i1,i2,i3)+&
!& (gmet(i1,1)*gmet(i2,2)*gmet(i3,1)+&
!& gmet(i1,2)*gmet(i2,1)*gmet(i3,1)+&
!& gmet(i1,1)*gmet(i2,1)*gmet(i3,2))*aa(ii,6)+&
!& (gmet(i1,1)*gmet(i2,2)*gmet(i3,2)+&
!& gmet(i1,2)*gmet(i2,1)*gmet(i3,2)+&
!& gmet(i1,2)*gmet(i2,2)*gmet(i3,1))*aa(ii,2)+&
!& (gmet(i1,1)*gmet(i2,3)*gmet(i3,1)+&
!& gmet(i1,3)*gmet(i2,1)*gmet(i3,1)+&
!& gmet(i1,1)*gmet(i2,1)*gmet(i3,3))*aa(ii,5)+&
!& (gmet(i1,1)*gmet(i2,3)*gmet(i3,3)+&
!& gmet(i1,3)*gmet(i2,1)*gmet(i3,3)+&
!& gmet(i1,3)*gmet(i2,3)*gmet(i3,1))*aa(ii,3)+&
!& (gmet(i1,2)*gmet(i2,2)*gmet(i3,3)+&
!& gmet(i1,2)*gmet(i2,3)*gmet(i3,2)+&
!& gmet(i1,3)*gmet(i2,2)*gmet(i3,2))*aa(ii,9)+&
!& (gmet(i1,2)*gmet(i2,3)*gmet(i3,3)+&
!& gmet(i1,3)*gmet(i2,2)*gmet(i3,3)+&
!& gmet(i1,3)*gmet(i2,3)*gmet(i3,2))*aa(ii,8)

!DEBUG
!write(std_out,*)' metcon : enter '
!stop
!ENDDEBUG
 if (rank==0) then
!  Simply copy scalar, re and im
   bb(1,1)=aa(1,1)
   bb(2,1)=aa(2,1)

 else if (rank==1) then
!  Apply gmet to input vector, re and im
   do ii=1,2
     do jj=1,3
       bb(ii,jj)=gmet(jj,1)*aa(ii,1)+gmet(jj,2)*aa(ii,2)+gmet(jj,3)*aa(ii,3)
     end do
   end do


 else if (rank==2) then
!  Apply rank2 expression, re and im
   do ii=1,2
!    Carry out g(c,d)*aa(c,d) contraction to get scalar
     scalar=gmet(1,1)*aa(ii,1)+gmet(2,2)*aa(ii,2)+&
&     gmet(3,3)*aa(ii,3)+2.0d0*(gmet(2,1)*aa(ii,6)+&
&     gmet(3,1)*aa(ii,5)+gmet(3,2)*aa(ii,4) )
!    Write out components of contraction
!    (1,1)->1
     bb(ii,1)=0.5d0*(3.0d0*(gmet(1,1)*gmet(1,1)*aa(ii,1)+&
&     gmet(1,2)*gmet(1,2)*aa(ii,2)+gmet(1,3)*gmet(1,3)*aa(ii,3)+&
&     (gmet(1,2)*gmet(1,1)+gmet(1,1)*gmet(1,2))*aa(ii,6)+&
&     (gmet(1,3)*gmet(1,1)+gmet(1,1)*gmet(1,3))*aa(ii,5)+&
&     (gmet(1,3)*gmet(1,2)+gmet(1,2)*gmet(1,3))*aa(ii,4) ) &
&     - gmet(1,1)*scalar)
!    (2,2)->2
     bb(ii,2)=0.5d0*(3.0d0*(gmet(2,1)*gmet(2,1)*aa(ii,1)+&
&     gmet(2,2)*gmet(2,2)*aa(ii,2)+gmet(2,3)*gmet(2,3)*aa(ii,3)+&
&     (gmet(2,2)*gmet(2,1)+gmet(2,1)*gmet(2,2))*aa(ii,6)+&
&     (gmet(2,3)*gmet(2,1)+gmet(2,1)*gmet(2,3))*aa(ii,5)+&
&     (gmet(2,3)*gmet(2,2)+gmet(2,2)*gmet(2,3))*aa(ii,4) )&
&     - gmet(2,2)*scalar)
!    (3,3)->3
     bb(ii,3)=0.5d0*(3.0d0*(gmet(3,1)*gmet(3,1)*aa(ii,1)+&
&     gmet(3,2)*gmet(3,2)*aa(ii,2)+gmet(3,3)*gmet(3,3)*aa(ii,3)+&
&     (gmet(3,2)*gmet(3,1)+gmet(3,1)*gmet(3,2))*aa(ii,6)+&
&     (gmet(3,3)*gmet(3,1)+gmet(3,1)*gmet(3,3))*aa(ii,5)+&
&     (gmet(3,3)*gmet(3,2)+gmet(3,2)*gmet(3,3))*aa(ii,4) )&
&     - gmet(3,3)*scalar)
!    (3,2)->4
     bb(ii,4)=0.5d0*(3.0d0*(gmet(3,1)*gmet(2,1)*aa(ii,1)+&
&     gmet(3,2)*gmet(2,2)*aa(ii,2)+gmet(3,3)*gmet(2,3)*aa(ii,3)+&
&     (gmet(3,2)*gmet(2,1)+gmet(3,1)*gmet(2,2))*aa(ii,6)+&
&     (gmet(3,3)*gmet(2,1)+gmet(3,1)*gmet(2,3))*aa(ii,5)+&
&     (gmet(3,3)*gmet(2,2)+gmet(3,2)*gmet(2,3))*aa(ii,4) )&
&     - gmet(3,2)*scalar)
!    (3,1)->5
     bb(ii,5)=0.5d0*(3.0d0*(gmet(3,1)*gmet(1,1)*aa(ii,1)+&
&     gmet(3,2)*gmet(1,2)*aa(ii,2)+gmet(3,3)*gmet(1,3)*aa(ii,3)+&
&     (gmet(3,2)*gmet(1,1)+gmet(3,1)*gmet(1,2))*aa(ii,6)+&
&     (gmet(3,3)*gmet(1,1)+gmet(3,1)*gmet(1,3))*aa(ii,5)+&
&     (gmet(3,3)*gmet(1,2)+gmet(3,2)*gmet(1,3))*aa(ii,4) )&
&     - gmet(3,1)*scalar)
!    (2,1)->6
     bb(ii,6)=0.5d0*(3.0d0*(gmet(2,1)*gmet(1,1)*aa(ii,1)+&
&     gmet(2,2)*gmet(1,2)*aa(ii,2)+gmet(2,3)*gmet(1,3)*aa(ii,3)+&
&     (gmet(2,2)*gmet(1,1)+gmet(2,1)*gmet(1,2))*aa(ii,6)+&
&     (gmet(2,3)*gmet(1,1)+gmet(2,1)*gmet(1,3))*aa(ii,5)+&
&     (gmet(2,3)*gmet(1,2)+gmet(2,2)*gmet(1,3))*aa(ii,4) ) &
&     - gmet(2,1)*scalar)
!    Include appropriate weights for multiplicity
     bb(ii,4)=2.d0*bb(ii,4)
     bb(ii,5)=2.d0*bb(ii,5)
     bb(ii,6)=2.d0*bb(ii,6)
   end do

 else if (rank==3) then
!  Apply rank2 expression, re and im
   do ii=1,2
!    Carry out g(l,m)g(j,n)*aa(l,m,n) contraction to get vector(j)
     do jj=1,3
       tmpiii=   gmet(1,1)*gmet(jj,1)*aa(ii,1)+&
&       gmet(2,2)*gmet(jj,2)*aa(ii,7)+&
&       gmet(3,3)*gmet(jj,3)*aa(ii,10)
       tmpijk=  (gmet(1,2)*gmet(jj,3)+&
&       gmet(3,1)*gmet(jj,2)+&
&       gmet(2,3)*gmet(jj,1)+&
&       gmet(3,2)*gmet(jj,1)+&
&       gmet(1,3)*gmet(jj,2)+&
&       gmet(2,1)*gmet(jj,3)) *aa(ii,4)
       vector(jj)=tmpiii + tmpijk +&
&       (gmet(1,2)*gmet(jj,1)+&
&       gmet(2,1)*gmet(jj,1)+&
&       gmet(1,1)*gmet(jj,2)) *aa(ii,6)+&
&       (gmet(1,2)*gmet(jj,2)+&
&       gmet(2,1)*gmet(jj,2)+&
&       gmet(2,2)*gmet(jj,1)) *aa(ii,2)+&
&       (gmet(1,3)*gmet(jj,1)+&
&       gmet(3,1)*gmet(jj,1)+&
&       gmet(1,1)*gmet(jj,3)) *aa(ii,5)+&
&       (gmet(1,3)*gmet(jj,3)+&
&       gmet(3,1)*gmet(jj,3)+&
&       gmet(3,3)*gmet(jj,1)) *aa(ii,3)+&
&       (gmet(2,3)*gmet(jj,2)+&
&       gmet(3,2)*gmet(jj,2)+&
&       gmet(2,2)*gmet(jj,3)) *aa(ii,9)+&
&       (gmet(2,3)*gmet(jj,3)+&
&       gmet(3,2)*gmet(jj,3)+&
&       gmet(3,3)*gmet(jj,2)) *aa(ii,8)
     end do
!    Write out components of contraction
!    (111)->1
     bb(ii,1) =2.5d0*con_met(ii,1,1,1)-1.5d0*(gmet(1,1)*vector(1))
!    (221)->2
     bb(ii,2) =2.5d0*con_met(ii,2,2,1)-0.5d0*(gmet(1,2)*vector(2)+&
&     gmet(1,2)*vector(2)+gmet(2,2)*vector(1))
!    (331)->3
     bb(ii,3) =2.5d0*con_met(ii,3,3,1)-0.5d0*(gmet(1,3)*vector(3)+&
&     gmet(1,3)*vector(3)+gmet(3,3)*vector(1))
!    (321)->4
     bb(ii,4) =2.5d0*con_met(ii,3,2,1)-0.5d0*(gmet(1,3)*vector(2)+&
&     gmet(1,2)*vector(3)+gmet(3,2)*vector(1))
!    (311)->5
     bb(ii,5) =2.5d0*con_met(ii,3,1,1)-0.5d0*(gmet(1,3)*vector(1)+&
&     gmet(1,1)*vector(3)+gmet(3,1)*vector(1))
!    (211)->6
     bb(ii,6) =2.5d0*con_met(ii,2,1,1)-0.5d0*(gmet(1,2)*vector(1)+&
&     gmet(1,1)*vector(2)+gmet(2,1)*vector(1))
!    (222)->7
     bb(ii,7) =2.5d0*con_met(ii,2,2,2)-1.5d0*(gmet(2,2)*vector(2))

!    (332)->8
     bb(ii,8) =2.5d0*con_met(ii,3,3,2)-0.5d0*(gmet(2,3)*vector(3)+&
&     gmet(2,3)*vector(3)+gmet(3,3)*vector(2))
!    (322)->9
     bb(ii,9) =2.5d0*con_met(ii,3,2,2)-0.5d0*(gmet(2,3)*vector(2)+&
&     gmet(2,2)*vector(3)+gmet(3,2)*vector(2))
!    (333)->10
     bb(ii,10)=2.5d0*con_met(ii,3,3,3)-1.5d0*(gmet(3,3)*vector(3))
!    Include appropriate weights for multiplicity
     bb(ii,2)=3.d0*bb(ii,2)
     bb(ii,3)=3.d0*bb(ii,3)
     bb(ii,4)=6.d0*bb(ii,4)
     bb(ii,5)=3.d0*bb(ii,5)
     bb(ii,6)=3.d0*bb(ii,6)
     bb(ii,8)=3.d0*bb(ii,8)
     bb(ii,9)=3.d0*bb(ii,9)
   end do

 else
   write(message, '(a,i0,a,a,a)' )&
&   'Input rank=',rank,' not allowed.',ch10,&
&   'Possible values are 0,1,2,3 only.'
   MSG_BUG(message)
 end if

 contains

   function con_met(ii,i1,i2,i3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'con_met'
!End of the abilint section

   real(dp) :: con_met
   integer :: ii,i1,i2,i3
   real(dp)::coniii,conijk

   coniii=gmet(i1,1)*gmet(i2,1)*gmet(i3,1)*aa(ii,1)+&
&   gmet(i1,2)*gmet(i2,2)*gmet(i3,2)*aa(ii,7)+&
&   gmet(i1,3)*gmet(i2,3)*gmet(i3,3)*aa(ii,10)
   conijk=aa(ii,4)*&
&   (gmet(i1,1)*gmet(i2,2)*gmet(i3,3)+&
&   gmet(i1,2)*gmet(i2,3)*gmet(i3,1)+&
&   gmet(i1,3)*gmet(i2,1)*gmet(i3,2)+&
&   gmet(i1,3)*gmet(i2,2)*gmet(i3,1)+&
&   gmet(i1,1)*gmet(i2,3)*gmet(i3,2)+&
&   gmet(i1,2)*gmet(i2,1)*gmet(i3,3))
   con_met=coniii+conijk+&
&   (gmet(i1,1)*gmet(i2,2)*gmet(i3,1)+&
&   gmet(i1,2)*gmet(i2,1)*gmet(i3,1)+&
&   gmet(i1,1)*gmet(i2,1)*gmet(i3,2))*aa(ii,6)+&
&   (gmet(i1,1)*gmet(i2,2)*gmet(i3,2)+&
&   gmet(i1,2)*gmet(i2,1)*gmet(i3,2)+&
&   gmet(i1,2)*gmet(i2,2)*gmet(i3,1))*aa(ii,2)+&
&   (gmet(i1,1)*gmet(i2,3)*gmet(i3,1)+&
&   gmet(i1,3)*gmet(i2,1)*gmet(i3,1)+&
&   gmet(i1,1)*gmet(i2,1)*gmet(i3,3))*aa(ii,5)+&
&   (gmet(i1,1)*gmet(i2,3)*gmet(i3,3)+&
&   gmet(i1,3)*gmet(i2,1)*gmet(i3,3)+&
&   gmet(i1,3)*gmet(i2,3)*gmet(i3,1))*aa(ii,3)+&
&   (gmet(i1,2)*gmet(i2,2)*gmet(i3,3)+&
&   gmet(i1,2)*gmet(i2,3)*gmet(i3,2)+&
&   gmet(i1,3)*gmet(i2,2)*gmet(i3,2))*aa(ii,9)+&
&   (gmet(i1,2)*gmet(i2,3)*gmet(i3,3)+&
&   gmet(i1,3)*gmet(i2,2)*gmet(i3,3)+&
&   gmet(i1,3)*gmet(i2,3)*gmet(i3,2))*aa(ii,8)

 end function con_met

end subroutine metcon
!!***

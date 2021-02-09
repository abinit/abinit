!!****m* ABINIT/m_metstr
!! NAME
!!  m_metstr
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group ()
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_metstr

 use defs_basis
 use m_errors
 use m_abicore

 implicit none

 private
!!***

 public :: metstr
!!***

contains
!!***

!!****f* ABINIT/metstr
!! NAME
!! metstr
!!
!! FUNCTION
!! Carries out specialized metric tensor operations needed for
!! the strain derivative of the l=0,1,2,3 nonlocal Kleinman-Bylander
!! pseudopotential operation.  Derivative is wrt a single (symmetric)
!! cartesian strain component.
!! Full advantage is taken of the full permutational symmetry of these
!! tensors.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (DRH, DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  istr=1,...6 specifies cartesian strain component 11,22,33,32,31,21
!!  rank=angular momentum
!!  iterm=1,2, or 3 as discussed below
!!  gmet(3,3)=metric tensor (array is symmetric but stored as 3x3)
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  aa(2,(rank+3)*(rank+4)/2)=unique elements of complex input tensor
!!    active size could be smaller, see rank discussion below
!!
!! OUTPUT
!!  bb(2,(rank+3)*(rank+4)/2)=unique elements of complex output tensor,
!!   rank discussed below
!!
!! NOTES
!! Based on metcon.f
!! All tensors are stored in a compressed storage mode defined below;
!! input and output conform to this scheme.
!! When tensor elements occur repeatedly due to symmetry, the
!! WEIGHT IS INCLUDED in the output tensor element to simplify later
!! contractions with other tensors of the same rank and form, i.e. the
!! next contraction is then simply a dot product over the unique elements.
!!
!! The matrix elements of the Kleinman-Bylander operator,
!!
!!      MKB(K\prim,K) = V_l(|K\prim|) P_l(K\prim,K) V_l(|K|)
!!
!! depend on the strain eps only through gmet.  Note that in the above
!! expession V_l=v_l/K^l, where v_l is the l component of the nonlocal
!! potential.  Also, P_l is a Legendre polynomial modified to be of
!! homogeneous order in K.  For example,
!!
!!      P_2 = (3(K\prim*K)^2 - |K\prim|^2|K|^2)/2
!!
!! Thus
!!
!!     (d/d eps)MKB(K\prim,K) = (d/d gmet)MKB(K\prim,K) (dgmet/d eps)
!!
!! has 3 terms coresponding to the 3 terms in the MKB product.  The rank
!! of the input(K) and output (K\prim) tensors are as follows for each term
!!
!! iterm=1  (d/d gmet)V_l(K\prim)     input rank = l,   output rank = l+2
!! iterm=2  (d/d gmet)P_l(K\prim,K)   input rank = l,   output rank = l
!! iterm=3  (d/d gmet)V_l(K)      input rank = l+2, output rank = l
!!
!! While playing a similar role to the routine metcon in caclulating
!! the coefficients to be used in constructing the output wavefunctions
!! of the Kleinman-Bylander operation, metstr can contract or expand its
!! tensor arguments, or neither.
!!
!! The compressed storage scheme is based on storing a symmetric 3x3 matrix as
!! $$
!!      \left( \begin{array}{ccc}
!!       1 & \cdots  & \cdots  \
!!       6 &    2    & \cdots  \
!!       5 &    4    &    3
!!       \end{array} \right)
!! $$
!!
!! PARENTS
!!      m_nonlop_pl
!!
!! CHILDREN
!!
!! SOURCE

subroutine metstr(istr,rank,iterm,gmet,gprimd,aa,bb)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istr,iterm,rank
!arrays
 real(dp),intent(in) :: aa(2,((rank+3)*(rank+4))/2),gmet(3,3),gprimd(3,3)
 real(dp),intent(out) :: bb(2,((rank+3)*(rank+4))/2)

!Local variables-------------------------------
!scalars
 integer,parameter :: mrank=3
 integer,save :: old_istr=0
 integer :: ii,jj,ka,kb,limitin,limitout,rankin,rankout
 character(len=500) :: message
!arrays
 integer,save :: cm_set(3,0:mrank),idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 real(dp),save :: cm(((mrank+3)*(mrank+4))/2,((mrank+3)*(mrank+4))/2,3,0:mrank)
 real(dp) :: dgmetds(3,3)

! *************************************************************************
 if (iterm <1 .or. iterm>3) then
   write(message, '(a,i0,a,a,a)' )&
&   'Input iterm=',iterm,' not allowed.',ch10,&
&   'Possible values are 1,2,3 only.'
   ABI_BUG(message)
 end if

 if(istr/=old_istr) then
   cm_set(:,:)=0
   old_istr=istr
 end if

 if(cm_set(iterm,rank)==0) then

   cm_set(iterm,rank)=1

   if(istr<1 .or. istr>6) then
     write(message,'(a,i0,a,a,a)')&
&     'Input istr=',istr,' not allowed.',ch10,&
&     'Possible values are 1,2,3,4,5,6 only.'
     ABI_BUG(message)
   end if

   ka=idx(2*istr-1);kb=idx(2*istr)
   do ii = 1,3
     dgmetds(:,ii)=-(gprimd(ka,:)*gprimd(kb,ii)+gprimd(kb,:)*gprimd(ka,ii))
   end do
!  For historical reasons:
   dgmetds(:,:)=0.5d0*dgmetds(:,:)

!
!  The code below was written by a Mathematica program and formatted by
!  a combination of editing scripts.  It is not intended to be read
!  by human beings, and certainly not to be modified by one.  Conceivably
!  it could be shortened somewhat by identifying common subexpressions.
!  However, it is only executed ONCE in each run for a given lattice and
!  strain component, so why worry.  Only a small double loop at the end
!  is executed on any but the first call.
!
   if (rank==0) then
     if(iterm==1) then
       cm(1,1,1,0)=dgmetds(1,1)
       cm(1,2,1,0)=dgmetds(2,2)
       cm(1,3,1,0)=dgmetds(3,3)
       cm(1,4,1,0)=2.d0*dgmetds(2,3)
       cm(1,5,1,0)=2.d0*dgmetds(1,3)
       cm(1,6,1,0)=2.d0*dgmetds(1,2)
     elseif(iterm==2) then
       cm(1,1,2,0)=0.d0
     elseif(iterm==3) then
       cm(1,1,3,0)= dgmetds(1,1)
       cm(2,1,3,0)= dgmetds(2,2)
       cm(3,1,3,0)= dgmetds(3,3)
       cm(4,1,3,0)= 2.d0*dgmetds(2,3)
       cm(5,1,3,0)= 2.d0*dgmetds(1,3)
       cm(6,1,3,0)= 2.d0*dgmetds(1,2)
     end if

   elseif(rank==1)then
     if(iterm==1)then
       cm(1,1,1,1)=gmet(1,1)*dgmetds(1,1)
       cm(2,1,1,1)=gmet(1,2)*dgmetds(1,1)
       cm(3,1,1,1)=gmet(1,3)*dgmetds(1,1)
       cm(1,2,1,1)=2*gmet(1,2)*dgmetds(1,2)+gmet(1,1)*dgmetds(2,2)
       cm(2,2,1,1)=2*gmet(2,2)*dgmetds(1,2)+gmet(1,2)*dgmetds(2,2)
       cm(3,2,1,1)=2*gmet(2,3)*dgmetds(1,2)+gmet(1,3)*dgmetds(2,2)
       cm(1,3,1,1)=2*gmet(1,3)*dgmetds(1,3)+gmet(1,1)*dgmetds(3,3)
       cm(2,3,1,1)=2*gmet(2,3)*dgmetds(1,3)+gmet(1,2)*dgmetds(3,3)
       cm(3,3,1,1)=2*gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3)
       cm(1,4,1,1)=2*(gmet(1,3)*dgmetds(1,2)+gmet(1,2)*dgmetds(1,3)+gmet(1,1)&
&       *dgmetds(2,3))
       cm(2,4,1,1)=2*(gmet(2,3)*dgmetds(1,2)+gmet(2,2)*dgmetds(1,3)+gmet(1,2)&
&       *dgmetds(2,3))
       cm(3,4,1,1)=2*(gmet(3,3)*dgmetds(1,2)+gmet(2,3)*dgmetds(1,3)+gmet(1,3)&
&       *dgmetds(2,3))
       cm(1,5,1,1)=gmet(1,3)*dgmetds(1,1)+2*gmet(1,1)*dgmetds(1,3)
       cm(2,5,1,1)=gmet(2,3)*dgmetds(1,1)+2*gmet(1,2)*dgmetds(1,3)
       cm(3,5,1,1)=gmet(3,3)*dgmetds(1,1)+2*gmet(1,3)*dgmetds(1,3)
       cm(1,6,1,1)=gmet(1,2)*dgmetds(1,1)+2*gmet(1,1)*dgmetds(1,2)
       cm(2,6,1,1)=gmet(2,2)*dgmetds(1,1)+2*gmet(1,2)*dgmetds(1,2)
       cm(3,6,1,1)=gmet(2,3)*dgmetds(1,1)+2*gmet(1,3)*dgmetds(1,2)
       cm(1,7,1,1)=gmet(1,2)*dgmetds(2,2)
       cm(2,7,1,1)=gmet(2,2)*dgmetds(2,2)
       cm(3,7,1,1)=gmet(2,3)*dgmetds(2,2)
       cm(1,8,1,1)=2*gmet(1,3)*dgmetds(2,3)+gmet(1,2)*dgmetds(3,3)
       cm(2,8,1,1)=2*gmet(2,3)*dgmetds(2,3)+gmet(2,2)*dgmetds(3,3)
       cm(3,8,1,1)=2*gmet(3,3)*dgmetds(2,3)+gmet(2,3)*dgmetds(3,3)
       cm(1,9,1,1)=gmet(1,3)*dgmetds(2,2)+2*gmet(1,2)*dgmetds(2,3)
       cm(2,9,1,1)=gmet(2,3)*dgmetds(2,2)+2*gmet(2,2)*dgmetds(2,3)
       cm(3,9,1,1)=gmet(3,3)*dgmetds(2,2)+2*gmet(2,3)*dgmetds(2,3)
       cm(1,10,1,1)=gmet(1,3)*dgmetds(3,3)
       cm(2,10,1,1)=gmet(2,3)*dgmetds(3,3)
       cm(3,10,1,1)=gmet(3,3)*dgmetds(3,3)
     elseif(iterm==2)then
       cm(1,1,2,1)=2*dgmetds(1,1)
       cm(2,1,2,1)=2*dgmetds(1,2)
       cm(3,1,2,1)=2*dgmetds(1,3)
       cm(1,2,2,1)=2*dgmetds(1,2)
       cm(2,2,2,1)=2*dgmetds(2,2)
       cm(3,2,2,1)=2*dgmetds(2,3)
       cm(1,3,2,1)=2*dgmetds(1,3)
       cm(2,3,2,1)=2*dgmetds(2,3)
       cm(3,3,2,1)=2*dgmetds(3,3)
     elseif(iterm==3)then
       cm(1,1,3,1)=gmet(1,1)*dgmetds(1,1)
       cm(2,1,3,1)=2*gmet(1,2)*dgmetds(1,2)+gmet(1,1)*dgmetds(2,2)
       cm(3,1,3,1)=2*gmet(1,3)*dgmetds(1,3)+gmet(1,1)*dgmetds(3,3)
       cm(4,1,3,1)=2*(gmet(1,3)*dgmetds(1,2)+gmet(1,2)*dgmetds(1,3)+gmet(1,1)&
&       *dgmetds(2,3))
       cm(5,1,3,1)=gmet(1,3)*dgmetds(1,1)+2*gmet(1,1)*dgmetds(1,3)
       cm(6,1,3,1)=gmet(1,2)*dgmetds(1,1)+2*gmet(1,1)*dgmetds(1,2)
       cm(7,1,3,1)=gmet(1,2)*dgmetds(2,2)
       cm(8,1,3,1)=2*gmet(1,3)*dgmetds(2,3)+gmet(1,2)*dgmetds(3,3)
       cm(9,1,3,1)=gmet(1,3)*dgmetds(2,2)+2*gmet(1,2)*dgmetds(2,3)
       cm(10,1,3,1)=gmet(1,3)*dgmetds(3,3)
       cm(1,2,3,1)=gmet(1,2)*dgmetds(1,1)
       cm(2,2,3,1)=2*gmet(2,2)*dgmetds(1,2)+gmet(1,2)*dgmetds(2,2)
       cm(3,2,3,1)=2*gmet(2,3)*dgmetds(1,3)+gmet(1,2)*dgmetds(3,3)
       cm(4,2,3,1)=2*(gmet(2,3)*dgmetds(1,2)+gmet(2,2)*dgmetds(1,3)+gmet(1,2)&
&       *dgmetds(2,3))
       cm(5,2,3,1)=gmet(2,3)*dgmetds(1,1)+2*gmet(1,2)*dgmetds(1,3)
       cm(6,2,3,1)=gmet(2,2)*dgmetds(1,1)+2*gmet(1,2)*dgmetds(1,2)
       cm(7,2,3,1)=gmet(2,2)*dgmetds(2,2)
       cm(8,2,3,1)=2*gmet(2,3)*dgmetds(2,3)+gmet(2,2)*dgmetds(3,3)
       cm(9,2,3,1)=gmet(2,3)*dgmetds(2,2)+2*gmet(2,2)*dgmetds(2,3)
       cm(10,2,3,1)=gmet(2,3)*dgmetds(3,3)
       cm(1,3,3,1)=gmet(1,3)*dgmetds(1,1)
       cm(2,3,3,1)=2*gmet(2,3)*dgmetds(1,2)+gmet(1,3)*dgmetds(2,2)
       cm(3,3,3,1)=2*gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3)
       cm(4,3,3,1)=2*(gmet(3,3)*dgmetds(1,2)+gmet(2,3)*dgmetds(1,3)+gmet(1,3)&
&       *dgmetds(2,3))
       cm(5,3,3,1)=gmet(3,3)*dgmetds(1,1)+2*gmet(1,3)*dgmetds(1,3)
       cm(6,3,3,1)=gmet(2,3)*dgmetds(1,1)+2*gmet(1,3)*dgmetds(1,2)
       cm(7,3,3,1)=gmet(2,3)*dgmetds(2,2)
       cm(8,3,3,1)=2*gmet(3,3)*dgmetds(2,3)+gmet(2,3)*dgmetds(3,3)
       cm(9,3,3,1)=gmet(3,3)*dgmetds(2,2)+2*gmet(2,3)*dgmetds(2,3)
       cm(10,3,3,1)=gmet(3,3)*dgmetds(3,3)
     end if

   elseif(rank==2)then
     if(iterm==1)then
       cm(1,1,1,2)=gmet(1,1)**2*dgmetds(1,1)
       cm(2,1,1,2)=((6*gmet(1,2)**2-2*gmet(1,1)*gmet(2,2))*dgmetds(1,1))&
&       /4.d0
       cm(3,1,1,2)=((6*gmet(1,3)**2-2*gmet(1,1)*gmet(3,3))*dgmetds(1,1))&
&       /4.d0
       cm(4,1,1,2)=((6*gmet(1,2)*gmet(1,3)-2*gmet(1,1)*gmet(2,3))*dgmetds(1,1))&
&       /2.d0
       cm(5,1,1,2)=2*gmet(1,1)*gmet(1,3)*dgmetds(1,1)
       cm(6,1,1,2)=2*gmet(1,1)*gmet(1,2)*dgmetds(1,1)
       cm(1,2,1,2)=1.5d0*gmet(1,2)**2*dgmetds(1,1)+4*gmet(1,1)*gmet(1,2)&
&       *dgmetds(1,2)+gmet(1,1)*(-0.5d0*gmet(2,2)*dgmetds(1,1)+gmet(1,1)&
&       *dgmetds(2,2))
       cm(2,2,1,2)=gmet(2,2)**2*dgmetds(1,1)+1.5d0*gmet(1,2)**2*dgmetds(2,2)&
&       +gmet(2,2)*(4*gmet(1,2)*dgmetds(1,2)-0.5d0*gmet(1,1)*dgmetds(2,2))
       cm(3,2,1,2)=1.5d0*gmet(2,3)**2*dgmetds(1,1)-0.5d0*gmet(2,2)*gmet(3,3)&
&       *dgmetds(1,1)+6*gmet(1,3)*gmet(2,3)*dgmetds(1,2)-2*gmet(1,2)&
&       *gmet(3,3)*dgmetds(1,2)+1.5d0*gmet(1,3)**2*dgmetds(2,2)-0.5d0*gmet(1,1)&
&       *gmet(3,3)*dgmetds(2,2)
       cm(4,2,1,2)=gmet(2,2)*(2*gmet(2,3)*dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,2))&
&       -gmet(1,1)*gmet(2,3)*dgmetds(2,2)+gmet(1,2)*(2*gmet(2,3)*dgmetds(1,2)&
&       +3*gmet(1,3)*dgmetds(2,2))
       cm(5,2,1,2)=gmet(2,3)*(3*gmet(1,2)*dgmetds(1,1)+6*gmet(1,1)*dgmetds(1,2))&
&       +gmet(1,3)*(-gmet(2,2)*dgmetds(1,1)+2*(gmet(1,2)*dgmetds(1,2)&
&       +gmet(1,1)*dgmetds(2,2)))
       cm(6,2,1,2)=2*gmet(1,2)**2*dgmetds(1,2)+6*gmet(1,1)*gmet(2,2)&
&       *dgmetds(1,2)+2*gmet(1,2)*(gmet(2,2)*dgmetds(1,1)+gmet(1,1)*dgmetds(2,2))
       cm(1,3,1,2)=1.5d0*gmet(1,3)**2*dgmetds(1,1)+4*gmet(1,1)*gmet(1,3)&
&       *dgmetds(1,3)+gmet(1,1)*(-0.5d0*gmet(3,3)*dgmetds(1,1)+gmet(1,1)&
&       *dgmetds(3,3))
       cm(2,3,1,2)=1.5d0*gmet(2,3)**2*dgmetds(1,1)+6*gmet(1,2)*gmet(2,3)&
&       *dgmetds(1,3)+1.5d0*gmet(1,2)**2*dgmetds(3,3)+gmet(2,2)*(-0.5d0*gmet(3,3)&
&       *dgmetds(1,1)-2*gmet(1,3)*dgmetds(1,3)-0.5d0*gmet(1,1)*dgmetds(3,3))
       cm(3,3,1,2)=gmet(3,3)**2*dgmetds(1,1)+1.5d0*gmet(1,3)**2*dgmetds(3,3)&
&       +gmet(3,3)*(4*gmet(1,3)*dgmetds(1,3)-0.5d0*gmet(1,1)*dgmetds(3,3))
       cm(4,3,1,2)=gmet(2,3)*(2*gmet(3,3)*dgmetds(1,1)+2*gmet(1,3)*dgmetds(1,3)&
&       -gmet(1,1)*dgmetds(3,3))+gmet(1,2)*(6*gmet(3,3)*dgmetds(1,3)&
&       +3*gmet(1,3)*dgmetds(3,3))
       cm(5,3,1,2)=2*gmet(1,3)**2*dgmetds(1,3)+6*gmet(1,1)*gmet(3,3)&
&       *dgmetds(1,3)+2*gmet(1,3)*(gmet(3,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(3,3))
       cm(6,3,1,2)=6*gmet(1,1)*gmet(2,3)*dgmetds(1,3)+gmet(1,3)*(3*gmet(2,3)&
&       *dgmetds(1,1)+2*gmet(1,2)*dgmetds(1,3))+gmet(1,2)*(-gmet(3,3)&
&       *dgmetds(1,1)+2*gmet(1,1)*dgmetds(3,3))
       cm(1,4,1,2)=gmet(1,2)*(3*gmet(1,3)*dgmetds(1,1)+4*gmet(1,1)*dgmetds(1,3))&
&       +gmet(1,1)*(-gmet(2,3)*dgmetds(1,1)+4*gmet(1,3)*dgmetds(1,2)&
&       +2*gmet(1,1)*dgmetds(2,3))
       cm(2,4,1,2)=gmet(2,2)*(2*gmet(2,3)*dgmetds(1,1)-2*gmet(1,3)*dgmetds(1,2)&
&       +4*gmet(1,2)*dgmetds(1,3)-gmet(1,1)*dgmetds(2,3))+gmet(1,2)*(6*gmet(2,3)&
&       *dgmetds(1,2)+3*gmet(1,2)*dgmetds(2,3))
       cm(3,4,1,2)=4*gmet(1,3)*gmet(3,3)*dgmetds(1,2)+gmet(2,3)*(2*gmet(3,3)&
&       *dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,3))+3*gmet(1,3)**2*dgmetds(2,3)&
&       +gmet(3,3)*(-2*gmet(1,2)*dgmetds(1,3)-gmet(1,1)*dgmetds(2,3))
       cm(4,4,1,2)=(2*gmet(2,3)**2*dgmetds(1,1)+gmet(2,2)*(6*gmet(3,3)&
&       *dgmetds(1,1)+12*gmet(1,3)*dgmetds(1,3))+gmet(2,3)*(4*gmet(1,3)&
&       *dgmetds(1,2)+4*gmet(1,2)*dgmetds(1,3)-4*gmet(1,1)*dgmetds(2,3))&
&       +12*gmet(1,2)*(gmet(3,3)*dgmetds(1,2)+gmet(1,3)*dgmetds(2,3)))&
&       /2.d0
       cm(5,4,1,2)=3*gmet(1,2)*gmet(3,3)*dgmetds(1,1)+2*gmet(1,3)**2*dgmetds(1,2)&
&       +6*gmet(1,1)*(gmet(3,3)*dgmetds(1,2)+gmet(2,3)*dgmetds(1,3))&
&       +gmet(1,3)*(1*gmet(2,3)*dgmetds(1,1)+2*gmet(1,2)*dgmetds(1,3)&
&       +4*gmet(1,1)*dgmetds(2,3))
       cm(6,4,1,2)=gmet(1,3)*(3*gmet(2,2)*dgmetds(1,1)+2*gmet(1,2)*dgmetds(1,2))&
&       +2*gmet(1,2)**2*dgmetds(1,3)+6*gmet(1,1)*(gmet(2,3)*dgmetds(1,2)&
&       +gmet(2,2)*dgmetds(1,3))+gmet(1,2)*(1*gmet(2,3)*dgmetds(1,1)&
&       +4*gmet(1,1)*dgmetds(2,3))
       cm(1,5,1,2)=2*gmet(1,1)*(gmet(1,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,3))
       cm(2,5,1,2)=-gmet(1,3)*gmet(2,2)*dgmetds(1,1)+3*gmet(1,2)*gmet(2,3)&
&       *dgmetds(1,1)+3*gmet(1,2)**2*dgmetds(1,3)-gmet(1,1)*gmet(2,2)&
&       *dgmetds(1,3)
       cm(3,5,1,2)=2*gmet(1,3)*gmet(3,3)*dgmetds(1,1)+3*gmet(1,3)**2*dgmetds(1,3)&
&       -gmet(1,1)*gmet(3,3)*dgmetds(1,3)
       cm(4,5,1,2)=3*gmet(1,2)*gmet(3,3)*dgmetds(1,1)-2*gmet(1,1)*gmet(2,3)&
&       *dgmetds(1,3)+gmet(1,3)*(1*gmet(2,3)*dgmetds(1,1)+6*gmet(1,2)&
&       *dgmetds(1,3))
       cm(5,5,1,2)=gmet(1,3)**2*dgmetds(1,1)+3*gmet(1,1)*gmet(3,3)*dgmetds(1,1)&
&       +4*gmet(1,1)*gmet(1,3)*dgmetds(1,3)
       cm(6,5,1,2)=3*gmet(1,1)*gmet(2,3)*dgmetds(1,1)+gmet(1,2)*(1*gmet(1,3)&
&       *dgmetds(1,1)+4*gmet(1,1)*dgmetds(1,3))
       cm(1,6,1,2)=2*gmet(1,1)*(gmet(1,2)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,2))
       cm(2,6,1,2)=2*gmet(1,2)*gmet(2,2)*dgmetds(1,1)+3*gmet(1,2)**2*dgmetds(1,2)&
&       -gmet(1,1)*gmet(2,2)*dgmetds(1,2)
       cm(3,6,1,2)=3*gmet(1,3)*gmet(2,3)*dgmetds(1,1)+3*gmet(1,3)**2*dgmetds(1,2)&
&       -gmet(3,3)*(gmet(1,2)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,2))
       cm(4,6,1,2)=gmet(2,3)*(1*gmet(1,2)*dgmetds(1,1)-2*gmet(1,1)*dgmetds(1,2))&
&       +gmet(1,3)*(3*gmet(2,2)*dgmetds(1,1)+6*gmet(1,2)*dgmetds(1,2))
       cm(5,6,1,2)=gmet(1,2)*gmet(1,3)*dgmetds(1,1)+gmet(1,1)*(3*gmet(2,3)&
&       *dgmetds(1,1)+4*gmet(1,3)*dgmetds(1,2))
       cm(6,6,1,2)=gmet(1,2)**2*dgmetds(1,1)+3*gmet(1,1)*gmet(2,2)*dgmetds(1,1)&
&       +4*gmet(1,1)*gmet(1,2)*dgmetds(1,2)
       cm(1,7,1,2)=3*gmet(1,2)**2*dgmetds(1,2)-gmet(1,1)*gmet(2,2)*dgmetds(1,2)&
&       +2*gmet(1,1)*gmet(1,2)*dgmetds(2,2)
       cm(2,7,1,2)=2*gmet(2,2)*(gmet(2,2)*dgmetds(1,2)+gmet(1,2)*dgmetds(2,2))
       cm(3,7,1,2)=3*gmet(2,3)**2*dgmetds(1,2)+3*gmet(1,3)*gmet(2,3)&
&       *dgmetds(2,2)-gmet(3,3)*(gmet(2,2)*dgmetds(1,2)+gmet(1,2)*dgmetds(2,2))
       cm(4,7,1,2)=gmet(1,2)*gmet(2,3)*dgmetds(2,2)+gmet(2,2)*(4*gmet(2,3)&
&       *dgmetds(1,2)+3*gmet(1,3)*dgmetds(2,2))
       cm(5,7,1,2)=gmet(2,3)*(6*gmet(1,2)*dgmetds(1,2)+3*gmet(1,1)*dgmetds(2,2))&
&       +gmet(1,3)*(-2*gmet(2,2)*dgmetds(1,2)+gmet(1,2)*dgmetds(2,2))
       cm(6,7,1,2)=4*gmet(1,2)*gmet(2,2)*dgmetds(1,2)+gmet(1,2)**2*dgmetds(2,2)&
&       +3*gmet(1,1)*gmet(2,2)*dgmetds(2,2)
       cm(1,8,1,2)=3*gmet(1,3)**2*dgmetds(1,2)+gmet(1,3)*(6*gmet(1,2)&
&       *dgmetds(1,3)+4*gmet(1,1)*dgmetds(2,3))+gmet(1,1)*(-gmet(3,3)&
&       *dgmetds(1,2)-2*gmet(2,3)*dgmetds(1,3)+2*gmet(1,2)*dgmetds(3,3))
       cm(2,8,1,2)=3*gmet(2,3)**2*dgmetds(1,2)+gmet(2,3)*(4*gmet(2,2)&
&       *dgmetds(1,3)+6*gmet(1,2)*dgmetds(2,3))+gmet(2,2)*(-gmet(3,3)&
&       *dgmetds(1,2)-2*gmet(1,3)*dgmetds(2,3)+2*gmet(1,2)*dgmetds(3,3))
       cm(3,8,1,2)=2*gmet(3,3)**2*dgmetds(1,2)+3*gmet(1,3)*gmet(2,3)&
&       *dgmetds(3,3)+gmet(3,3)*(4*gmet(2,3)*dgmetds(1,3)+4*gmet(1,3)&
&       *dgmetds(2,3)-gmet(1,2)*dgmetds(3,3))
       cm(4,8,1,2)=2*gmet(2,3)**2*dgmetds(1,3)+6*gmet(1,2)*gmet(3,3)&
&       *dgmetds(2,3)+gmet(2,3)*(4*gmet(3,3)*dgmetds(1,2)+2*gmet(1,3)&
&       *dgmetds(2,3)+gmet(1,2)*dgmetds(3,3))+gmet(2,2)*(6*gmet(3,3)&
&       *dgmetds(1,3)+3*gmet(1,3)*dgmetds(3,3))
       cm(5,8,1,2)=6*gmet(1,2)*gmet(3,3)*dgmetds(1,3)+2*gmet(1,3)**2*dgmetds(2,3)&
&       +gmet(1,3)*(4*gmet(3,3)*dgmetds(1,2)+2*gmet(2,3)*dgmetds(1,3)&
&       +gmet(1,2)*dgmetds(3,3))+gmet(1,1)*(6*gmet(3,3)*dgmetds(2,3)&
&       +3*gmet(2,3)*dgmetds(3,3))
       cm(6,8,1,2)=(2*(6*gmet(1,3)*gmet(2,3)-2*gmet(1,2)*gmet(3,3))*dgmetds(1,2)&
&       +4*(3*gmet(1,3)*gmet(2,2)+gmet(1,2)*gmet(2,3))*dgmetds(1,3)+4*(1*gmet(1,2)&
&       *gmet(1,3)+3*gmet(1,1)*gmet(2,3))*dgmetds(2,3)+2*(1*gmet(1,2)&
&       **2+3*gmet(1,1)*gmet(2,2))*dgmetds(3,3))/2.d0
       cm(1,9,1,2)=3*gmet(1,2)**2*dgmetds(1,3)+gmet(1,1)*(-2*gmet(2,3)&
&       *dgmetds(1,2)-gmet(2,2)*dgmetds(1,3)+2*gmet(1,3)*dgmetds(2,2))&
&       +gmet(1,2)*(6*gmet(1,3)*dgmetds(1,2)+4*gmet(1,1)*dgmetds(2,3))
       cm(2,9,1,2)=2*gmet(2,2)**2*dgmetds(1,3)+3*gmet(1,2)*gmet(2,3)&
&       *dgmetds(2,2)+gmet(2,2)*(4*gmet(2,3)*dgmetds(1,2)-gmet(1,3)*dgmetds(2,2)&
&       +4*gmet(1,2)*dgmetds(2,3))
       cm(3,9,1,2)=3*gmet(2,3)**2*dgmetds(1,3)+gmet(3,3)*(-gmet(2,2)&
&       *dgmetds(1,3)+2*gmet(1,3)*dgmetds(2,2)-2*gmet(1,2)*dgmetds(2,3))&
&       +gmet(2,3)*(4*gmet(3,3)*dgmetds(1,2)+6*gmet(1,3)*dgmetds(2,3))
       cm(4,9,1,2)=2*gmet(2,3)**2*dgmetds(1,2)+3*gmet(1,2)*gmet(3,3)&
&       *dgmetds(2,2)+gmet(2,3)*(4*gmet(2,2)*dgmetds(1,3)+gmet(1,3)*dgmetds(2,2)&
&       +2*gmet(1,2)*dgmetds(2,3))+6*gmet(2,2)*(gmet(3,3)*dgmetds(1,2)&
&       +gmet(1,3)*dgmetds(2,3))
       cm(5,9,1,2)=(12*gmet(1,2)*(gmet(3,3)*dgmetds(1,2)+gmet(2,3)*dgmetds(1,3))&
&       +2*gmet(1,3)**2*dgmetds(2,2)+gmet(1,3)*(4*gmet(2,3)*dgmetds(1,2)&
&       -4*gmet(2,2)*dgmetds(1,3)+4*gmet(1,2)*dgmetds(2,3))+gmet(1,1)&
&       *(6*gmet(3,3)*dgmetds(2,2)+12*gmet(2,3)*dgmetds(2,3)))/2.d0
       cm(6,9,1,2)=gmet(1,2)*(2*gmet(2,3)*dgmetds(1,2)+4*gmet(2,2)*dgmetds(1,3))&
&       +gmet(1,3)*(6*gmet(2,2)*dgmetds(1,2)+gmet(1,2)*dgmetds(2,2))&
&       +2*gmet(1,2)**2*dgmetds(2,3)+gmet(1,1)*(3*gmet(2,3)*dgmetds(2,2)&
&       +6*gmet(2,2)*dgmetds(2,3))
       cm(1,10,1,2)=3*gmet(1,3)**2*dgmetds(1,3)-gmet(1,1)*gmet(3,3)*dgmetds(1,3)&
&       +2*gmet(1,1)*gmet(1,3)*dgmetds(3,3)
       cm(2,10,1,2)=3*gmet(2,3)**2*dgmetds(1,3)+3*gmet(1,2)*gmet(2,3)&
&       *dgmetds(3,3)-gmet(2,2)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(3,10,1,2)=2*gmet(3,3)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(4,10,1,2)=3*gmet(1,2)*gmet(3,3)*dgmetds(3,3)+gmet(2,3)*(4*gmet(3,3)&
&       *dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(5,10,1,2)=4*gmet(1,3)*gmet(3,3)*dgmetds(1,3)+gmet(1,3)**2*dgmetds(3,3)&
&       +3*gmet(1,1)*gmet(3,3)*dgmetds(3,3)
       cm(6,10,1,2)=-2*gmet(1,2)*gmet(3,3)*dgmetds(1,3)+3*gmet(1,1)*gmet(2,3)&
&       *dgmetds(3,3)+gmet(1,3)*(6*gmet(2,3)*dgmetds(1,3)+gmet(1,2)*dgmetds(3,3))
       cm(1,11,1,2)=((6*gmet(1,2)**2-2*gmet(1,1)*gmet(2,2))*dgmetds(2,2))&
&       /4.d0
       cm(2,11,1,2)=gmet(2,2)**2*dgmetds(2,2)
       cm(3,11,1,2)=((6*gmet(2,3)**2-2*gmet(2,2)*gmet(3,3))*dgmetds(2,2))&
&       /4.d0
       cm(4,11,1,2)=2*gmet(2,2)*gmet(2,3)*dgmetds(2,2)
       cm(5,11,1,2)=((-2*gmet(1,3)*gmet(2,2)+6*gmet(1,2)*gmet(2,3))*dgmetds(2,2))&
&       /2.d0
       cm(6,11,1,2)=2*gmet(1,2)*gmet(2,2)*dgmetds(2,2)
       cm(1,12,1,2)=1.5d0*gmet(1,3)**2*dgmetds(2,2)+6*gmet(1,2)*gmet(1,3)&
&       *dgmetds(2,3)+1.5d0*gmet(1,2)**2*dgmetds(3,3)+gmet(1,1)*(-0.5d0*gmet(3,3)&
&       *dgmetds(2,2)-2*gmet(2,3)*dgmetds(2,3)-0.5d0*gmet(2,2)*dgmetds(3,3))
       cm(2,12,1,2)=1.5d0*gmet(2,3)**2*dgmetds(2,2)+4*gmet(2,2)*gmet(2,3)&
&       *dgmetds(2,3)+gmet(2,2)*(-0.5d0*gmet(3,3)*dgmetds(2,2)+gmet(2,2)&
&       *dgmetds(3,3))
       cm(3,12,1,2)=gmet(3,3)**2*dgmetds(2,2)+1.5d0*gmet(2,3)**2*dgmetds(3,3)&
&       +gmet(3,3)*(4*gmet(2,3)*dgmetds(2,3)-0.5d0*gmet(2,2)*dgmetds(3,3))
       cm(4,12,1,2)=2*gmet(2,3)**2*dgmetds(2,3)+6*gmet(2,2)*gmet(3,3)&
&       *dgmetds(2,3)+2*gmet(2,3)*(gmet(3,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(3,3))
       cm(5,12,1,2)=gmet(1,3)*(2*gmet(3,3)*dgmetds(2,2)+2*gmet(2,3)*dgmetds(2,3)&
&       -gmet(2,2)*dgmetds(3,3))+gmet(1,2)*(6*gmet(3,3)*dgmetds(2,3)&
&       +3*gmet(2,3)*dgmetds(3,3))
       cm(6,12,1,2)=gmet(1,3)*(3*gmet(2,3)*dgmetds(2,2)+6*gmet(2,2)*dgmetds(2,3))&
&       +gmet(1,2)*(-gmet(3,3)*dgmetds(2,2)+2*(gmet(2,3)*dgmetds(2,3)&
&       +gmet(2,2)*dgmetds(3,3)))
       cm(1,13,1,2)=3*gmet(1,2)*gmet(1,3)*dgmetds(2,2)+3*gmet(1,2)**2*dgmetds(2,3)&
&       -gmet(1,1)*(gmet(2,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(2,3))
       cm(2,13,1,2)=2*gmet(2,2)*(gmet(2,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(2,3))
       cm(3,13,1,2)=2*gmet(2,3)*gmet(3,3)*dgmetds(2,2)+3*gmet(2,3)**2*dgmetds(2,3)&
&       -gmet(2,2)*gmet(3,3)*dgmetds(2,3)
       cm(4,13,1,2)=gmet(2,3)**2*dgmetds(2,2)+3*gmet(2,2)*gmet(3,3)*dgmetds(2,2)&
&       +4*gmet(2,2)*gmet(2,3)*dgmetds(2,3)
       cm(5,13,1,2)=gmet(1,3)*(1*gmet(2,3)*dgmetds(2,2)-2*gmet(2,2)*dgmetds(2,3))&
&       +gmet(1,2)*(3*gmet(3,3)*dgmetds(2,2)+6*gmet(2,3)*dgmetds(2,3))
       cm(6,13,1,2)=3*gmet(1,3)*gmet(2,2)*dgmetds(2,2)+gmet(1,2)*(1*gmet(2,3)&
&       *dgmetds(2,2)+4*gmet(2,2)*dgmetds(2,3))
       cm(1,14,1,2)=3*gmet(1,3)**2*dgmetds(2,3)+3*gmet(1,2)*gmet(1,3)&
&       *dgmetds(3,3)-gmet(1,1)*(gmet(3,3)*dgmetds(2,3)+gmet(2,3)*dgmetds(3,3))
       cm(2,14,1,2)=3*gmet(2,3)**2*dgmetds(2,3)-gmet(2,2)*gmet(3,3)*dgmetds(2,3)&
&       +2*gmet(2,2)*gmet(2,3)*dgmetds(3,3)
       cm(3,14,1,2)=2*gmet(3,3)*(gmet(3,3)*dgmetds(2,3)+gmet(2,3)*dgmetds(3,3))
       cm(4,14,1,2)=4*gmet(2,3)*gmet(3,3)*dgmetds(2,3)+gmet(2,3)**2*dgmetds(3,3)&
&       +3*gmet(2,2)*gmet(3,3)*dgmetds(3,3)
       cm(5,14,1,2)=3*gmet(1,2)*gmet(3,3)*dgmetds(3,3)+gmet(1,3)*(4*gmet(3,3)&
&       *dgmetds(2,3)+gmet(2,3)*dgmetds(3,3))
       cm(6,14,1,2)=gmet(1,3)*(6*gmet(2,3)*dgmetds(2,3)+3*gmet(2,2)*dgmetds(3,3))&
&       +gmet(1,2)*(-2*gmet(3,3)*dgmetds(2,3)+gmet(2,3)*dgmetds(3,3))
       cm(1,15,1,2)=((6*gmet(1,3)**2-2*gmet(1,1)*gmet(3,3))*dgmetds(3,3))&
&       /4.d0
       cm(2,15,1,2)=((6*gmet(2,3)**2-2*gmet(2,2)*gmet(3,3))*dgmetds(3,3))&
&       /4.d0
       cm(3,15,1,2)=gmet(3,3)**2*dgmetds(3,3)
       cm(4,15,1,2)=2*gmet(2,3)*gmet(3,3)*dgmetds(3,3)
       cm(5,15,1,2)=2*gmet(1,3)*gmet(3,3)*dgmetds(3,3)
       cm(6,15,1,2)=((6*gmet(1,3)*gmet(2,3)-2*gmet(1,2)*gmet(3,3))*dgmetds(3,3))&
&       /2.d0
     elseif(iterm==2)then
       cm(1,1,2,2)=4*gmet(1,1)*dgmetds(1,1)
       cm(2,1,2,2)=-gmet(2,2)*dgmetds(1,1)+6*gmet(1,2)*dgmetds(1,2)-gmet(1,1)&
&       *dgmetds(2,2)
       cm(3,1,2,2)=-gmet(3,3)*dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,3)-gmet(1,1)&
&       *dgmetds(3,3)
       cm(4,1,2,2)=-2*gmet(2,3)*dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,2)&
&       +6*gmet(1,2)*dgmetds(1,3)-2*gmet(1,1)*dgmetds(2,3)
       cm(5,1,2,2)=4*(gmet(1,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,3))
       cm(6,1,2,2)=4*(gmet(1,2)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,2))
       cm(1,2,2,2)=-gmet(2,2)*dgmetds(1,1)+6*gmet(1,2)*dgmetds(1,2)-gmet(1,1)&
&       *dgmetds(2,2)
       cm(2,2,2,2)=4*gmet(2,2)*dgmetds(2,2)
       cm(3,2,2,2)=-gmet(3,3)*dgmetds(2,2)+6*gmet(2,3)*dgmetds(2,3)-gmet(2,2)&
&       *dgmetds(3,3)
       cm(4,2,2,2)=4*(gmet(2,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(2,3))
       cm(5,2,2,2)=6*gmet(2,3)*dgmetds(1,2)-2*gmet(2,2)*dgmetds(1,3)&
&       -2*gmet(1,3)*dgmetds(2,2)+6*gmet(1,2)*dgmetds(2,3)
       cm(6,2,2,2)=4*(gmet(2,2)*dgmetds(1,2)+gmet(1,2)*dgmetds(2,2))
       cm(1,3,2,2)=-gmet(3,3)*dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,3)-gmet(1,1)&
&       *dgmetds(3,3)
       cm(2,3,2,2)=-gmet(3,3)*dgmetds(2,2)+6*gmet(2,3)*dgmetds(2,3)-gmet(2,2)&
&       *dgmetds(3,3)
       cm(3,3,2,2)=4*gmet(3,3)*dgmetds(3,3)
       cm(4,3,2,2)=4*(gmet(3,3)*dgmetds(2,3)+gmet(2,3)*dgmetds(3,3))
       cm(5,3,2,2)=4*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(6,3,2,2)=-2*gmet(3,3)*dgmetds(1,2)+6*gmet(2,3)*dgmetds(1,3)&
&       +6*gmet(1,3)*dgmetds(2,3)-2*gmet(1,2)*dgmetds(3,3)
       cm(1,4,2,2)=-2*gmet(2,3)*dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,2)&
&       +6*gmet(1,2)*dgmetds(1,3)-2*gmet(1,1)*dgmetds(2,3)
       cm(2,4,2,2)=4*(gmet(2,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(2,3))
       cm(3,4,2,2)=4*(gmet(3,3)*dgmetds(2,3)+gmet(2,3)*dgmetds(3,3))
       cm(4,4,2,2)=6*gmet(3,3)*dgmetds(2,2)+4*gmet(2,3)*dgmetds(2,3)&
&       +6*gmet(2,2)*dgmetds(3,3)
       cm(5,4,2,2)=6*gmet(3,3)*dgmetds(1,2)+2*gmet(2,3)*dgmetds(1,3)&
&       +2*gmet(1,3)*dgmetds(2,3)+6*gmet(1,2)*dgmetds(3,3)
       cm(6,4,2,2)=2*gmet(2,3)*dgmetds(1,2)+6*gmet(2,2)*dgmetds(1,3)&
&       +6*gmet(1,3)*dgmetds(2,2)+2*gmet(1,2)*dgmetds(2,3)
       cm(1,5,2,2)=4*(gmet(1,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,3))
       cm(2,5,2,2)=6*gmet(2,3)*dgmetds(1,2)-2*gmet(2,2)*dgmetds(1,3)&
&       -2*gmet(1,3)*dgmetds(2,2)+6*gmet(1,2)*dgmetds(2,3)
       cm(3,5,2,2)=4*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(4,5,2,2)=6*gmet(3,3)*dgmetds(1,2)+2*gmet(2,3)*dgmetds(1,3)&
&       +2*gmet(1,3)*dgmetds(2,3)+6*gmet(1,2)*dgmetds(3,3)
       cm(5,5,2,2)=6*gmet(3,3)*dgmetds(1,1)+4*gmet(1,3)*dgmetds(1,3)&
&       +6*gmet(1,1)*dgmetds(3,3)
       cm(6,5,2,2)=6*gmet(2,3)*dgmetds(1,1)+2*gmet(1,3)*dgmetds(1,2)&
&       +2*gmet(1,2)*dgmetds(1,3)+6*gmet(1,1)*dgmetds(2,3)
       cm(1,6,2,2)=4*(gmet(1,2)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,2))
       cm(2,6,2,2)=4*(gmet(2,2)*dgmetds(1,2)+gmet(1,2)*dgmetds(2,2))
       cm(3,6,2,2)=-2*gmet(3,3)*dgmetds(1,2)+6*gmet(2,3)*dgmetds(1,3)&
&       +6*gmet(1,3)*dgmetds(2,3)-2*gmet(1,2)*dgmetds(3,3)
       cm(4,6,2,2)=2*gmet(2,3)*dgmetds(1,2)+6*gmet(2,2)*dgmetds(1,3)&
&       +6*gmet(1,3)*dgmetds(2,2)+2*gmet(1,2)*dgmetds(2,3)
       cm(5,6,2,2)=6*gmet(2,3)*dgmetds(1,1)+2*gmet(1,3)*dgmetds(1,2)&
&       +2*gmet(1,2)*dgmetds(1,3)+6*gmet(1,1)*dgmetds(2,3)
       cm(6,6,2,2)=6*gmet(2,2)*dgmetds(1,1)+4*gmet(1,2)*dgmetds(1,2)&
&       +6*gmet(1,1)*dgmetds(2,2)
     elseif(iterm==3)then
       cm(1,1,3,2)=gmet(1,1)**2*dgmetds(1,1)
       cm(2,1,3,2)=1.5d0*gmet(1,2)**2*dgmetds(1,1)+4*gmet(1,1)*gmet(1,2)&
&       *dgmetds(1,2)+gmet(1,1)*(-0.5d0*gmet(2,2)*dgmetds(1,1)+gmet(1,1)&
&       *dgmetds(2,2))
       cm(3,1,3,2)=1.5d0*gmet(1,3)**2*dgmetds(1,1)+4*gmet(1,1)*gmet(1,3)&
&       *dgmetds(1,3)+gmet(1,1)*(-0.5d0*gmet(3,3)*dgmetds(1,1)+gmet(1,1)&
&       *dgmetds(3,3))
       cm(4,1,3,2)=gmet(1,2)*(3*gmet(1,3)*dgmetds(1,1)+4*gmet(1,1)*dgmetds(1,3))&
&       +gmet(1,1)*(-gmet(2,3)*dgmetds(1,1)+4*gmet(1,3)*dgmetds(1,2)&
&       +2*gmet(1,1)*dgmetds(2,3))
       cm(5,1,3,2)=2*gmet(1,1)*(gmet(1,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,3))
       cm(6,1,3,2)=2*gmet(1,1)*(gmet(1,2)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,2))
       cm(7,1,3,2)=3*gmet(1,2)**2*dgmetds(1,2)-gmet(1,1)*gmet(2,2)*dgmetds(1,2)&
&       +2*gmet(1,1)*gmet(1,2)*dgmetds(2,2)
       cm(8,1,3,2)=3*gmet(1,3)**2*dgmetds(1,2)+gmet(1,3)*(6*gmet(1,2)&
&       *dgmetds(1,3)+4*gmet(1,1)*dgmetds(2,3))+gmet(1,1)*(-gmet(3,3)&
&       *dgmetds(1,2)-2*gmet(2,3)*dgmetds(1,3)+2*gmet(1,2)*dgmetds(3,3))
       cm(9,1,3,2)=3*gmet(1,2)**2*dgmetds(1,3)+gmet(1,1)*(-2*gmet(2,3)&
&       *dgmetds(1,2)-gmet(2,2)*dgmetds(1,3)+2*gmet(1,3)*dgmetds(2,2))&
&       +gmet(1,2)*(6*gmet(1,3)*dgmetds(1,2)+4*gmet(1,1)*dgmetds(2,3))
       cm(10,1,3,2)=3*gmet(1,3)**2*dgmetds(1,3)-gmet(1,1)*gmet(3,3)*dgmetds(1,3)&
&       +2*gmet(1,1)*gmet(1,3)*dgmetds(3,3)
       cm(11,1,3,2)=((6*gmet(1,2)**2-2*gmet(1,1)*gmet(2,2))*dgmetds(2,2))&
&       /4.d0
       cm(12,1,3,2)=1.5d0*gmet(1,3)**2*dgmetds(2,2)+6*gmet(1,2)*gmet(1,3)&
&       *dgmetds(2,3)+1.5d0*gmet(1,2)**2*dgmetds(3,3)+gmet(1,1)*(-0.5d0*gmet(3,3)&
&       *dgmetds(2,2)-2*gmet(2,3)*dgmetds(2,3)-0.5d0*gmet(2,2)*dgmetds(3,3))
       cm(13,1,3,2)=3*gmet(1,2)*gmet(1,3)*dgmetds(2,2)+3*gmet(1,2)**2*dgmetds(2,3)&
&       -gmet(1,1)*(gmet(2,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(2,3))
       cm(14,1,3,2)=3*gmet(1,3)**2*dgmetds(2,3)+3*gmet(1,2)*gmet(1,3)&
&       *dgmetds(3,3)-gmet(1,1)*(gmet(3,3)*dgmetds(2,3)+gmet(2,3)*dgmetds(3,3))
       cm(15,1,3,2)=((6*gmet(1,3)**2-2*gmet(1,1)*gmet(3,3))*dgmetds(3,3))&
&       /4.d0
       cm(1,2,3,2)=((6*gmet(1,2)**2-2*gmet(1,1)*gmet(2,2))*dgmetds(1,1))&
&       /4.d0
       cm(2,2,3,2)=gmet(2,2)**2*dgmetds(1,1)+1.5d0*gmet(1,2)**2*dgmetds(2,2)&
&       +gmet(2,2)*(4*gmet(1,2)*dgmetds(1,2)-0.5d0*gmet(1,1)*dgmetds(2,2))
       cm(3,2,3,2)=1.5d0*gmet(2,3)**2*dgmetds(1,1)+6*gmet(1,2)*gmet(2,3)&
&       *dgmetds(1,3)+1.5d0*gmet(1,2)**2*dgmetds(3,3)+gmet(2,2)*(-0.5d0*gmet(3,3)&
&       *dgmetds(1,1)-2*gmet(1,3)*dgmetds(1,3)-0.5d0*gmet(1,1)*dgmetds(3,3))
       cm(4,2,3,2)=gmet(2,2)*(2*gmet(2,3)*dgmetds(1,1)-2*gmet(1,3)*dgmetds(1,2)&
&       +4*gmet(1,2)*dgmetds(1,3)-gmet(1,1)*dgmetds(2,3))+gmet(1,2)*(6*gmet(2,3)&
&       *dgmetds(1,2)+3*gmet(1,2)*dgmetds(2,3))
       cm(5,2,3,2)=-gmet(1,3)*gmet(2,2)*dgmetds(1,1)+3*gmet(1,2)*gmet(2,3)&
&       *dgmetds(1,1)+3*gmet(1,2)**2*dgmetds(1,3)-gmet(1,1)*gmet(2,2)&
&       *dgmetds(1,3)
       cm(6,2,3,2)=2*gmet(1,2)*gmet(2,2)*dgmetds(1,1)+3*gmet(1,2)**2*dgmetds(1,2)&
&       -gmet(1,1)*gmet(2,2)*dgmetds(1,2)
       cm(7,2,3,2)=2*gmet(2,2)*(gmet(2,2)*dgmetds(1,2)+gmet(1,2)*dgmetds(2,2))
       cm(8,2,3,2)=3*gmet(2,3)**2*dgmetds(1,2)+gmet(2,3)*(4*gmet(2,2)&
&       *dgmetds(1,3)+6*gmet(1,2)*dgmetds(2,3))+gmet(2,2)*(-gmet(3,3)&
&       *dgmetds(1,2)-2*gmet(1,3)*dgmetds(2,3)+2*gmet(1,2)*dgmetds(3,3))
       cm(9,2,3,2)=2*gmet(2,2)**2*dgmetds(1,3)+3*gmet(1,2)*gmet(2,3)&
&       *dgmetds(2,2)+gmet(2,2)*(4*gmet(2,3)*dgmetds(1,2)-gmet(1,3)*dgmetds(2,2)&
&       +4*gmet(1,2)*dgmetds(2,3))
       cm(10,2,3,2)=3*gmet(2,3)**2*dgmetds(1,3)+3*gmet(1,2)*gmet(2,3)&
&       *dgmetds(3,3)-gmet(2,2)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(11,2,3,2)=gmet(2,2)**2*dgmetds(2,2)
       cm(12,2,3,2)=1.5d0*gmet(2,3)**2*dgmetds(2,2)+4*gmet(2,2)*gmet(2,3)&
&       *dgmetds(2,3)+gmet(2,2)*(-0.5d0*gmet(3,3)*dgmetds(2,2)+gmet(2,2)&
&       *dgmetds(3,3))
       cm(13,2,3,2)=2*gmet(2,2)*(gmet(2,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(2,3))
       cm(14,2,3,2)=3*gmet(2,3)**2*dgmetds(2,3)-gmet(2,2)*gmet(3,3)*dgmetds(2,3)&
&       +2*gmet(2,2)*gmet(2,3)*dgmetds(3,3)
       cm(15,2,3,2)=((6*gmet(2,3)**2-2*gmet(2,2)*gmet(3,3))*dgmetds(3,3))&
&       /4.d0
       cm(1,3,3,2)=((6*gmet(1,3)**2-2*gmet(1,1)*gmet(3,3))*dgmetds(1,1))&
&       /4.d0
       cm(2,3,3,2)=1.5d0*gmet(2,3)**2*dgmetds(1,1)-0.5d0*gmet(2,2)*gmet(3,3)&
&       *dgmetds(1,1)+6*gmet(1,3)*gmet(2,3)*dgmetds(1,2)-2*gmet(1,2)&
&       *gmet(3,3)*dgmetds(1,2)+1.5d0*gmet(1,3)**2*dgmetds(2,2)-0.5d0*gmet(1,1)&
&       *gmet(3,3)*dgmetds(2,2)
       cm(3,3,3,2)=gmet(3,3)**2*dgmetds(1,1)+1.5d0*gmet(1,3)**2*dgmetds(3,3)&
&       +gmet(3,3)*(4*gmet(1,3)*dgmetds(1,3)-0.5d0*gmet(1,1)*dgmetds(3,3))
       cm(4,3,3,2)=4*gmet(1,3)*gmet(3,3)*dgmetds(1,2)+gmet(2,3)*(2*gmet(3,3)&
&       *dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,3))+3*gmet(1,3)**2*dgmetds(2,3)&
&       +gmet(3,3)*(-2*gmet(1,2)*dgmetds(1,3)-gmet(1,1)*dgmetds(2,3))
       cm(5,3,3,2)=2*gmet(1,3)*gmet(3,3)*dgmetds(1,1)+3*gmet(1,3)**2*dgmetds(1,3)&
&       -gmet(1,1)*gmet(3,3)*dgmetds(1,3)
       cm(6,3,3,2)=3*gmet(1,3)*gmet(2,3)*dgmetds(1,1)+3*gmet(1,3)**2*dgmetds(1,2)&
&       -gmet(3,3)*(gmet(1,2)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,2))
       cm(7,3,3,2)=3*gmet(2,3)**2*dgmetds(1,2)+3*gmet(1,3)*gmet(2,3)&
&       *dgmetds(2,2)-gmet(3,3)*(gmet(2,2)*dgmetds(1,2)+gmet(1,2)*dgmetds(2,2))
       cm(8,3,3,2)=2*gmet(3,3)**2*dgmetds(1,2)+3*gmet(1,3)*gmet(2,3)&
&       *dgmetds(3,3)+gmet(3,3)*(4*gmet(2,3)*dgmetds(1,3)+4*gmet(1,3)&
&       *dgmetds(2,3)-gmet(1,2)*dgmetds(3,3))
       cm(9,3,3,2)=3*gmet(2,3)**2*dgmetds(1,3)+gmet(3,3)*(-gmet(2,2)&
&       *dgmetds(1,3)+2*gmet(1,3)*dgmetds(2,2)-2*gmet(1,2)*dgmetds(2,3))&
&       +gmet(2,3)*(4*gmet(3,3)*dgmetds(1,2)+6*gmet(1,3)*dgmetds(2,3))
       cm(10,3,3,2)=2*gmet(3,3)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(11,3,3,2)=((6*gmet(2,3)**2-2*gmet(2,2)*gmet(3,3))*dgmetds(2,2))&
&       /4.d0
       cm(12,3,3,2)=gmet(3,3)**2*dgmetds(2,2)+1.5d0*gmet(2,3)**2*dgmetds(3,3)&
&       +gmet(3,3)*(4*gmet(2,3)*dgmetds(2,3)-0.5d0*gmet(2,2)*dgmetds(3,3))
       cm(13,3,3,2)=2*gmet(2,3)*gmet(3,3)*dgmetds(2,2)+3*gmet(2,3)**2*dgmetds(2,3)&
&       -gmet(2,2)*gmet(3,3)*dgmetds(2,3)
       cm(14,3,3,2)=2*gmet(3,3)*(gmet(3,3)*dgmetds(2,3)+gmet(2,3)*dgmetds(3,3))
       cm(15,3,3,2)=gmet(3,3)**2*dgmetds(3,3)
       cm(1,4,3,2)=((6*gmet(1,2)*gmet(1,3)-2*gmet(1,1)*gmet(2,3))*dgmetds(1,1))&
&       /2.d0
       cm(2,4,3,2)=gmet(2,2)*(2*gmet(2,3)*dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,2))&
&       -gmet(1,1)*gmet(2,3)*dgmetds(2,2)+gmet(1,2)*(2*gmet(2,3)*dgmetds(1,2)&
&       +3*gmet(1,3)*dgmetds(2,2))
       cm(3,4,3,2)=gmet(2,3)*(2*gmet(3,3)*dgmetds(1,1)+2*gmet(1,3)*dgmetds(1,3)&
&       -gmet(1,1)*dgmetds(3,3))+gmet(1,2)*(6*gmet(3,3)*dgmetds(1,3)&
&       +3*gmet(1,3)*dgmetds(3,3))
       cm(4,4,3,2)=(2*gmet(2,3)**2*dgmetds(1,1)+gmet(2,2)*(6*gmet(3,3)&
&       *dgmetds(1,1)+12*gmet(1,3)*dgmetds(1,3))+gmet(2,3)*(4*gmet(1,3)&
&       *dgmetds(1,2)+4*gmet(1,2)*dgmetds(1,3)-4*gmet(1,1)*dgmetds(2,3))&
&       +12*gmet(1,2)*(gmet(3,3)*dgmetds(1,2)+gmet(1,3)*dgmetds(2,3)))&
&       /2.d0
       cm(5,4,3,2)=3*gmet(1,2)*gmet(3,3)*dgmetds(1,1)-2*gmet(1,1)*gmet(2,3)&
&       *dgmetds(1,3)+gmet(1,3)*(1*gmet(2,3)*dgmetds(1,1)+6*gmet(1,2)&
&       *dgmetds(1,3))
       cm(6,4,3,2)=gmet(2,3)*(1*gmet(1,2)*dgmetds(1,1)-2*gmet(1,1)*dgmetds(1,2))&
&       +gmet(1,3)*(3*gmet(2,2)*dgmetds(1,1)+6*gmet(1,2)*dgmetds(1,2))
       cm(7,4,3,2)=gmet(1,2)*gmet(2,3)*dgmetds(2,2)+gmet(2,2)*(4*gmet(2,3)&
&       *dgmetds(1,2)+3*gmet(1,3)*dgmetds(2,2))
       cm(8,4,3,2)=2*gmet(2,3)**2*dgmetds(1,3)+6*gmet(1,2)*gmet(3,3)&
&       *dgmetds(2,3)+gmet(2,3)*(4*gmet(3,3)*dgmetds(1,2)+2*gmet(1,3)&
&       *dgmetds(2,3)+gmet(1,2)*dgmetds(3,3))+gmet(2,2)*(6*gmet(3,3)&
&       *dgmetds(1,3)+3*gmet(1,3)*dgmetds(3,3))
       cm(9,4,3,2)=2*gmet(2,3)**2*dgmetds(1,2)+3*gmet(1,2)*gmet(3,3)&
&       *dgmetds(2,2)+gmet(2,3)*(4*gmet(2,2)*dgmetds(1,3)+gmet(1,3)*dgmetds(2,2)&
&       +2*gmet(1,2)*dgmetds(2,3))+6*gmet(2,2)*(gmet(3,3)*dgmetds(1,2)&
&       +gmet(1,3)*dgmetds(2,3))
       cm(10,4,3,2)=3*gmet(1,2)*gmet(3,3)*dgmetds(3,3)+gmet(2,3)*(4*gmet(3,3)&
&       *dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(11,4,3,2)=2*gmet(2,2)*gmet(2,3)*dgmetds(2,2)
       cm(12,4,3,2)=2*gmet(2,3)**2*dgmetds(2,3)+6*gmet(2,2)*gmet(3,3)&
&       *dgmetds(2,3)+2*gmet(2,3)*(gmet(3,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(3,3))
       cm(13,4,3,2)=gmet(2,3)**2*dgmetds(2,2)+3*gmet(2,2)*gmet(3,3)*dgmetds(2,2)&
&       +4*gmet(2,2)*gmet(2,3)*dgmetds(2,3)
       cm(14,4,3,2)=4*gmet(2,3)*gmet(3,3)*dgmetds(2,3)+gmet(2,3)**2*dgmetds(3,3)&
&       +3*gmet(2,2)*gmet(3,3)*dgmetds(3,3)
       cm(15,4,3,2)=2*gmet(2,3)*gmet(3,3)*dgmetds(3,3)
       cm(1,5,3,2)=2*gmet(1,1)*gmet(1,3)*dgmetds(1,1)
       cm(2,5,3,2)=gmet(2,3)*(3*gmet(1,2)*dgmetds(1,1)+6*gmet(1,1)*dgmetds(1,2))&
&       +gmet(1,3)*(-gmet(2,2)*dgmetds(1,1)+2*(gmet(1,2)*dgmetds(1,2)&
&       +gmet(1,1)*dgmetds(2,2)))
       cm(3,5,3,2)=2*gmet(1,3)**2*dgmetds(1,3)+6*gmet(1,1)*gmet(3,3)&
&       *dgmetds(1,3)+2*gmet(1,3)*(gmet(3,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(3,3))
       cm(4,5,3,2)=3*gmet(1,2)*gmet(3,3)*dgmetds(1,1)+2*gmet(1,3)**2*dgmetds(1,2)&
&       +6*gmet(1,1)*(gmet(3,3)*dgmetds(1,2)+gmet(2,3)*dgmetds(1,3))&
&       +gmet(1,3)*(1*gmet(2,3)*dgmetds(1,1)+2*gmet(1,2)*dgmetds(1,3)&
&       +4*gmet(1,1)*dgmetds(2,3))
       cm(5,5,3,2)=gmet(1,3)**2*dgmetds(1,1)+3*gmet(1,1)*gmet(3,3)*dgmetds(1,1)&
&       +4*gmet(1,1)*gmet(1,3)*dgmetds(1,3)
       cm(6,5,3,2)=gmet(1,2)*gmet(1,3)*dgmetds(1,1)+gmet(1,1)*(3*gmet(2,3)&
&       *dgmetds(1,1)+4*gmet(1,3)*dgmetds(1,2))
       cm(7,5,3,2)=gmet(2,3)*(6*gmet(1,2)*dgmetds(1,2)+3*gmet(1,1)*dgmetds(2,2))&
&       +gmet(1,3)*(-2*gmet(2,2)*dgmetds(1,2)+gmet(1,2)*dgmetds(2,2))
       cm(8,5,3,2)=6*gmet(1,2)*gmet(3,3)*dgmetds(1,3)+2*gmet(1,3)**2*dgmetds(2,3)&
&       +gmet(1,3)*(4*gmet(3,3)*dgmetds(1,2)+2*gmet(2,3)*dgmetds(1,3)&
&       +gmet(1,2)*dgmetds(3,3))+gmet(1,1)*(6*gmet(3,3)*dgmetds(2,3)&
&       +3*gmet(2,3)*dgmetds(3,3))
       cm(9,5,3,2)=(12*gmet(1,2)*(gmet(3,3)*dgmetds(1,2)+gmet(2,3)*dgmetds(1,3))&
&       +2*gmet(1,3)**2*dgmetds(2,2)+gmet(1,3)*(4*gmet(2,3)*dgmetds(1,2)&
&       -4*gmet(2,2)*dgmetds(1,3)+4*gmet(1,2)*dgmetds(2,3))+gmet(1,1)&
&       *(6*gmet(3,3)*dgmetds(2,2)+12*gmet(2,3)*dgmetds(2,3)))/2.d0
       cm(10,5,3,2)=4*gmet(1,3)*gmet(3,3)*dgmetds(1,3)+gmet(1,3)**2*dgmetds(3,3)&
&       +3*gmet(1,1)*gmet(3,3)*dgmetds(3,3)
       cm(11,5,3,2)=((-2*gmet(1,3)*gmet(2,2)+6*gmet(1,2)*gmet(2,3))*dgmetds(2,2))&
&       /2.d0
       cm(12,5,3,2)=gmet(1,3)*(2*gmet(3,3)*dgmetds(2,2)+2*gmet(2,3)*dgmetds(2,3)&
&       -gmet(2,2)*dgmetds(3,3))+gmet(1,2)*(6*gmet(3,3)*dgmetds(2,3)&
&       +3*gmet(2,3)*dgmetds(3,3))
       cm(13,5,3,2)=gmet(1,3)*(1*gmet(2,3)*dgmetds(2,2)-2*gmet(2,2)*dgmetds(2,3))&
&       +gmet(1,2)*(3*gmet(3,3)*dgmetds(2,2)+6*gmet(2,3)*dgmetds(2,3))
       cm(14,5,3,2)=3*gmet(1,2)*gmet(3,3)*dgmetds(3,3)+gmet(1,3)*(4*gmet(3,3)&
&       *dgmetds(2,3)+gmet(2,3)*dgmetds(3,3))
       cm(15,5,3,2)=2*gmet(1,3)*gmet(3,3)*dgmetds(3,3)
       cm(1,6,3,2)=2*gmet(1,1)*gmet(1,2)*dgmetds(1,1)
       cm(2,6,3,2)=2*gmet(1,2)**2*dgmetds(1,2)+6*gmet(1,1)*gmet(2,2)&
&       *dgmetds(1,2)+2*gmet(1,2)*(gmet(2,2)*dgmetds(1,1)+gmet(1,1)*dgmetds(2,2))
       cm(3,6,3,2)=6*gmet(1,1)*gmet(2,3)*dgmetds(1,3)+gmet(1,3)*(3*gmet(2,3)&
&       *dgmetds(1,1)+2*gmet(1,2)*dgmetds(1,3))+gmet(1,2)*(-gmet(3,3)&
&       *dgmetds(1,1)+2*gmet(1,1)*dgmetds(3,3))
       cm(4,6,3,2)=gmet(1,3)*(3*gmet(2,2)*dgmetds(1,1)+2*gmet(1,2)*dgmetds(1,2))&
&       +2*gmet(1,2)**2*dgmetds(1,3)+6*gmet(1,1)*(gmet(2,3)*dgmetds(1,2)&
&       +gmet(2,2)*dgmetds(1,3))+gmet(1,2)*(1*gmet(2,3)*dgmetds(1,1)&
&       +4*gmet(1,1)*dgmetds(2,3))
       cm(5,6,3,2)=3*gmet(1,1)*gmet(2,3)*dgmetds(1,1)+gmet(1,2)*(1*gmet(1,3)&
&       *dgmetds(1,1)+4*gmet(1,1)*dgmetds(1,3))
       cm(6,6,3,2)=gmet(1,2)**2*dgmetds(1,1)+3*gmet(1,1)*gmet(2,2)*dgmetds(1,1)&
&       +4*gmet(1,1)*gmet(1,2)*dgmetds(1,2)
       cm(7,6,3,2)=4*gmet(1,2)*gmet(2,2)*dgmetds(1,2)+gmet(1,2)**2*dgmetds(2,2)&
&       +3*gmet(1,1)*gmet(2,2)*dgmetds(2,2)
       cm(8,6,3,2)=(2*(6*gmet(1,3)*gmet(2,3)-2*gmet(1,2)*gmet(3,3))*dgmetds(1,2)&
&       +4*(3*gmet(1,3)*gmet(2,2)+gmet(1,2)*gmet(2,3))*dgmetds(1,3)+4*(1*gmet(1,2)&
&       *gmet(1,3)+3*gmet(1,1)*gmet(2,3))*dgmetds(2,3)+2*(1*gmet(1,2)&
&       **2+3*gmet(1,1)*gmet(2,2))*dgmetds(3,3))/2.d0
       cm(9,6,3,2)=gmet(1,2)*(2*gmet(2,3)*dgmetds(1,2)+4*gmet(2,2)*dgmetds(1,3))&
&       +gmet(1,3)*(6*gmet(2,2)*dgmetds(1,2)+gmet(1,2)*dgmetds(2,2))&
&       +2*gmet(1,2)**2*dgmetds(2,3)+gmet(1,1)*(3*gmet(2,3)*dgmetds(2,2)&
&       +6*gmet(2,2)*dgmetds(2,3))
       cm(10,6,3,2)=-2*gmet(1,2)*gmet(3,3)*dgmetds(1,3)+3*gmet(1,1)*gmet(2,3)&
&       *dgmetds(3,3)+gmet(1,3)*(6*gmet(2,3)*dgmetds(1,3)+gmet(1,2)*dgmetds(3,3))
       cm(11,6,3,2)=2*gmet(1,2)*gmet(2,2)*dgmetds(2,2)
       cm(12,6,3,2)=gmet(1,3)*(3*gmet(2,3)*dgmetds(2,2)+6*gmet(2,2)*dgmetds(2,3))&
&       +gmet(1,2)*(-gmet(3,3)*dgmetds(2,2)+2*(gmet(2,3)*dgmetds(2,3)&
&       +gmet(2,2)*dgmetds(3,3)))
       cm(13,6,3,2)=3*gmet(1,3)*gmet(2,2)*dgmetds(2,2)+gmet(1,2)*(1*gmet(2,3)&
&       *dgmetds(2,2)+4*gmet(2,2)*dgmetds(2,3))
       cm(14,6,3,2)=gmet(1,3)*(6*gmet(2,3)*dgmetds(2,3)+3*gmet(2,2)*dgmetds(3,3))&
&       +gmet(1,2)*(-2*gmet(3,3)*dgmetds(2,3)+gmet(2,3)*dgmetds(3,3))
       cm(15,6,3,2)=((6*gmet(1,3)*gmet(2,3)-2*gmet(1,2)*gmet(3,3))*dgmetds(3,3))&
&       /2.d0
     end if

   elseif(rank==3)then
     if(iterm==1)then
       cm(1,1,1,3)=gmet(1,1)**3*dgmetds(1,1)
       cm(2,1,1,3)=gmet(1,1)*(4.5d0*gmet(1,2)**2-1.5d0*gmet(1,1)*gmet(2,2))&
&       *dgmetds(1,1)
       cm(3,1,1,3)=gmet(1,1)*(4.5d0*gmet(1,3)**2-1.5d0*gmet(1,1)*gmet(3,3))&
&       *dgmetds(1,1)
       cm(4,1,1,3)=(gmet(1,1)*(54*gmet(1,2)*gmet(1,3)-18*gmet(1,1)*gmet(2,3))&
&       *dgmetds(1,1))/6.d0
       cm(5,1,1,3)=3*gmet(1,1)**2*gmet(1,3)*dgmetds(1,1)
       cm(6,1,1,3)=3*gmet(1,1)**2*gmet(1,2)*dgmetds(1,1)
       cm(7,1,1,3)=gmet(1,2)*(2.5d0*gmet(1,2)**2-1.5d0*gmet(1,1)*gmet(2,2))&
&       *dgmetds(1,1)
       cm(8,1,1,3)=((-36*gmet(1,1)*gmet(1,3)*gmet(2,3)+gmet(1,2)*(90*gmet(1,3)&
&       **2-18*gmet(1,1)*gmet(3,3)))*dgmetds(1,1))/12.d0
       cm(9,1,1,3)=((90*gmet(1,2)**2*gmet(1,3)-18*gmet(1,1)*gmet(1,3)&
&       *gmet(2,2)-36*gmet(1,1)*gmet(1,2)*gmet(2,3))*dgmetds(1,1))/12.d0
       cm(10,1,1,3)=gmet(1,3)*(2.5d0*gmet(1,3)**2-1.5d0*gmet(1,1)*gmet(3,3))&
&       *dgmetds(1,1)
       cm(1,2,1,3)=gmet(1,1)*(4.5d0*gmet(1,2)**2*dgmetds(1,1)+6*gmet(1,1)&
&       *gmet(1,2)*dgmetds(1,2)+gmet(1,1)*(-1.5d0*gmet(2,2)*dgmetds(1,1)&
&       +gmet(1,1)*dgmetds(2,2)))
       cm(2,2,1,3)=3*gmet(1,2)**3*dgmetds(1,2)+15*gmet(1,1)*gmet(1,2)&
&       *gmet(2,2)*dgmetds(1,2)+gmet(1,1)*gmet(2,2)*(6*gmet(2,2)*dgmetds(1,1)&
&       -1.5d0*gmet(1,1)*dgmetds(2,2))+gmet(1,2)**2*(3*gmet(2,2)*dgmetds(1,1)&
&       +4.5d0*gmet(1,1)*dgmetds(2,2))
       cm(3,2,1,3)=-3*gmet(1,2)**2*gmet(3,3)*dgmetds(1,1)+gmet(1,3)*gmet(2,3)&
&       *(9*gmet(1,2)*dgmetds(1,1)+24*gmet(1,1)*dgmetds(1,2))+gmet(1,1)&
&       *(7.5d0*gmet(2,3)**2*dgmetds(1,1)-1.5d0*gmet(2,2)*gmet(3,3)*dgmetds(1,1)&
&       -9*gmet(1,2)*gmet(3,3)*dgmetds(1,2))-1.5d0*gmet(1,1)**2*gmet(3,3)&
&       *dgmetds(2,2)+gmet(1,3)**2*(-3*gmet(2,2)*dgmetds(1,1)+3*gmet(1,2)&
&       *dgmetds(1,2)+4.5d0*gmet(1,1)*dgmetds(2,2))
       cm(4,2,1,3)=gmet(1,2)**2*(3*gmet(2,3)*dgmetds(1,1)+6*gmet(1,3)&
&       *dgmetds(1,2))+gmet(1,1)*(gmet(2,2)*(12*gmet(2,3)*dgmetds(1,1)&
&       +24*gmet(1,3)*dgmetds(1,2))-3*gmet(1,1)*gmet(2,3)*dgmetds(2,2))&
&       +gmet(1,2)*(6*gmet(1,1)*gmet(2,3)*dgmetds(1,2)+gmet(1,3)*(3*gmet(2,2)&
&       *dgmetds(1,1)+9*gmet(1,1)*dgmetds(2,2)))
       cm(5,2,1,3)=1.5d0*gmet(1,2)**2*gmet(1,3)*dgmetds(1,1)+gmet(1,1)&
&       *gmet(1,2)*(12*gmet(2,3)*dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,2))&
&       +gmet(1,1)*(12*gmet(1,1)*gmet(2,3)*dgmetds(1,2)+gmet(1,3)*(-4.5d0*gmet(2,2)&
&       *dgmetds(1,1)+3*gmet(1,1)*dgmetds(2,2)))
       cm(6,2,1,3)=1.5d0*gmet(1,2)**3*dgmetds(1,1)+6*gmet(1,1)*gmet(1,2)&
&       **2*dgmetds(1,2)+12*gmet(1,1)**2*gmet(2,2)*dgmetds(1,2)+gmet(1,1)&
&       *gmet(1,2)*(7.5d0*gmet(2,2)*dgmetds(1,1)+3*gmet(1,1)*dgmetds(2,2))
       cm(7,2,1,3)=9*gmet(1,2)**2*gmet(2,2)*dgmetds(1,2)-3*gmet(1,1)&
&       *gmet(2,2)**2*dgmetds(1,2)+2.5d0*gmet(1,2)**3*dgmetds(2,2)+gmet(1,2)&
&       *gmet(2,2)*(3*gmet(2,2)*dgmetds(1,1)-1.5d0*gmet(1,1)*dgmetds(2,2))
       cm(8,2,1,3)=(6*(48*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(6*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+12*(30*gmet(1,3)**2*gmet(2,2)&
&       +36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)&
&       *(-12*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+2*(-36*gmet(1,1)&
&       *gmet(1,3)*gmet(2,3)+gmet(1,2)*(90*gmet(1,3)**2-18*gmet(1,1)&
&       *gmet(3,3)))*dgmetds(2,2))/24.d0
       cm(9,2,1,3)=gmet(2,3)*(3*gmet(1,2)**2*dgmetds(1,2)-9*gmet(1,1)&
&       *gmet(2,2)*dgmetds(1,2)+gmet(1,2)*(3*gmet(2,2)*dgmetds(1,1)-3*gmet(1,1)&
&       *dgmetds(2,2)))+gmet(1,3)*(6*gmet(2,2)**2*dgmetds(1,1)+7.5d0*gmet(1,2)&
&       **2*dgmetds(2,2)+gmet(2,2)*(24*gmet(1,2)*dgmetds(1,2)-1.5d0*gmet(1,1)&
&       *dgmetds(2,2)))
       cm(10,2,1,3)=(1080*gmet(1,3)**2*gmet(2,3)*dgmetds(1,2)-216*gmet(2,3)&
&       *gmet(3,3)*(gmet(1,2)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,2))+180*gmet(1,3)&
&       **3*dgmetds(2,2)+gmet(1,3)*(540*gmet(2,3)**2*dgmetds(1,1)+gmet(3,3)&
&       *(-108*gmet(2,2)*dgmetds(1,1)-432*gmet(1,2)*dgmetds(1,2)-108*gmet(1,1)&
&       *dgmetds(2,2))))/72.d0
       cm(1,3,1,3)=gmet(1,1)*(4.5d0*gmet(1,3)**2*dgmetds(1,1)+6*gmet(1,1)&
&       *gmet(1,3)*dgmetds(1,3)+gmet(1,1)*(-1.5d0*gmet(3,3)*dgmetds(1,1)&
&       +gmet(1,1)*dgmetds(3,3)))
       cm(2,3,1,3)=-3*gmet(1,3)**2*gmet(2,2)*dgmetds(1,1)-3*gmet(1,2)&
&       **2*gmet(3,3)*dgmetds(1,1)+gmet(1,3)*(9*gmet(1,2)*gmet(2,3)*dgmetds(1,1)&
&       +3*gmet(1,2)**2*dgmetds(1,3)-9*gmet(1,1)*gmet(2,2)*dgmetds(1,3))&
&       -1.5d0*gmet(1,1)**2*gmet(2,2)*dgmetds(3,3)+gmet(1,1)*(7.5d0*gmet(2,3)&
&       **2*dgmetds(1,1)-1.5d0*gmet(2,2)*gmet(3,3)*dgmetds(1,1)+24*gmet(1,2)&
&       *gmet(2,3)*dgmetds(1,3)+4.5d0*gmet(1,2)**2*dgmetds(3,3))
       cm(3,3,1,3)=3*gmet(1,3)**3*dgmetds(1,3)+15*gmet(1,1)*gmet(1,3)&
&       *gmet(3,3)*dgmetds(1,3)+gmet(1,1)*gmet(3,3)*(6*gmet(3,3)*dgmetds(1,1)&
&       -1.5d0*gmet(1,1)*dgmetds(3,3))+gmet(1,3)**2*(3*gmet(3,3)*dgmetds(1,1)&
&       +4.5d0*gmet(1,1)*dgmetds(3,3))
       cm(4,3,1,3)=gmet(1,3)**2*(3*gmet(2,3)*dgmetds(1,1)+6*gmet(1,2)&
&       *dgmetds(1,3))+gmet(1,1)*(24*gmet(1,2)*gmet(3,3)*dgmetds(1,3)&
&       +gmet(2,3)*(12*gmet(3,3)*dgmetds(1,1)-3*gmet(1,1)*dgmetds(3,3)))&
&       +gmet(1,3)*(6*gmet(1,1)*gmet(2,3)*dgmetds(1,3)+gmet(1,2)*(3*gmet(3,3)&
&       *dgmetds(1,1)+9*gmet(1,1)*dgmetds(3,3)))
       cm(5,3,1,3)=1.5d0*gmet(1,3)**3*dgmetds(1,1)+6*gmet(1,1)*gmet(1,3)&
&       **2*dgmetds(1,3)+12*gmet(1,1)**2*gmet(3,3)*dgmetds(1,3)+gmet(1,1)&
&       *gmet(1,3)*(7.5d0*gmet(3,3)*dgmetds(1,1)+3*gmet(1,1)*dgmetds(3,3))
       cm(6,3,1,3)=12*gmet(1,1)*gmet(2,3)*(gmet(1,3)*dgmetds(1,1)+gmet(1,1)&
&       *dgmetds(1,3))+gmet(1,2)*(1.5d0*gmet(1,3)**2*dgmetds(1,1)+6*gmet(1,1)&
&       *gmet(1,3)*dgmetds(1,3)+gmet(1,1)*(-4.5d0*gmet(3,3)*dgmetds(1,1)&
&       +3*gmet(1,1)*dgmetds(3,3)))
       cm(7,3,1,3)=(6*(-36*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+12*(-36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,2)+90*gmet(1,2)**2*gmet(2,3)-18*gmet(1,1)*gmet(2,2)*gmet(2,3))&
&       *dgmetds(1,3)+2*(90*gmet(1,2)**3-54*gmet(1,1)*gmet(1,2)*gmet(2,2))&
&       *dgmetds(3,3))/72.d0
       cm(8,3,1,3)=gmet(1,3)**2*(3*gmet(2,3)*dgmetds(1,3)+7.5d0*gmet(1,2)&
&       *dgmetds(3,3))+gmet(1,3)*(24*gmet(1,2)*gmet(3,3)*dgmetds(1,3)&
&       +gmet(2,3)*(3*gmet(3,3)*dgmetds(1,1)-3*gmet(1,1)*dgmetds(3,3)))&
&       +gmet(3,3)*(-9*gmet(1,1)*gmet(2,3)*dgmetds(1,3)+gmet(1,2)*(6*gmet(3,3)&
&       *dgmetds(1,1)-1.5d0*gmet(1,1)*dgmetds(3,3)))
       cm(9,3,1,3)=(6*(48*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(6*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+12*(-12*gmet(1,3)**2*gmet(2,2)&
&       +36*gmet(1,2)*gmet(1,3)*gmet(2,3)+30*gmet(1,2)**2*gmet(3,3)+gmet(1,1)&
&       *(-12*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,3)+2*(90*gmet(1,2)&
&       **2*gmet(1,3)-18*gmet(1,1)*gmet(1,3)*gmet(2,2)-36*gmet(1,1)*gmet(1,2)&
&       *gmet(2,3))*dgmetds(3,3))/24.d0
       cm(10,3,1,3)=9*gmet(1,3)**2*gmet(3,3)*dgmetds(1,3)-3*gmet(1,1)&
&       *gmet(3,3)**2*dgmetds(1,3)+2.5d0*gmet(1,3)**3*dgmetds(3,3)+gmet(1,3)&
&       *gmet(3,3)*(3*gmet(3,3)*dgmetds(1,1)-1.5d0*gmet(1,1)*dgmetds(3,3))
       cm(1,4,1,3)=gmet(1,1)*(gmet(1,2)*(9*gmet(1,3)*dgmetds(1,1)+6*gmet(1,1)&
&       *dgmetds(1,3))+gmet(1,1)*(-3*gmet(2,3)*dgmetds(1,1)+6*gmet(1,3)&
&       *dgmetds(1,2)+2*gmet(1,1)*dgmetds(2,3)))
       cm(2,4,1,3)=3*gmet(1,2)**3*dgmetds(1,3)+gmet(1,2)*(3*gmet(1,3)&
&       *gmet(2,2)*dgmetds(1,1)+gmet(1,1)*(24*gmet(2,3)*dgmetds(1,2)&
&       +15*gmet(2,2)*dgmetds(1,3)))+gmet(1,1)*gmet(2,2)*(12*gmet(2,3)&
&       *dgmetds(1,1)-9*gmet(1,3)*dgmetds(1,2)-3*gmet(1,1)*dgmetds(2,3))&
&       +gmet(1,2)**2*(3*gmet(2,3)*dgmetds(1,1)+3*gmet(1,3)*dgmetds(1,2)&
&       +9*gmet(1,1)*dgmetds(2,3))
       cm(3,4,1,3)=3*gmet(1,3)**3*dgmetds(1,2)+gmet(1,3)*(3*gmet(1,2)&
&       *gmet(3,3)*dgmetds(1,1)+gmet(1,1)*(15*gmet(3,3)*dgmetds(1,2)&
&       +24*gmet(2,3)*dgmetds(1,3)))+gmet(1,1)*gmet(3,3)*(12*gmet(2,3)&
&       *dgmetds(1,1)-9*gmet(1,2)*dgmetds(1,3)-3*gmet(1,1)*dgmetds(2,3))&
&       +gmet(1,3)**2*(3*gmet(2,3)*dgmetds(1,1)+3*gmet(1,2)*dgmetds(1,3)&
&       +9*gmet(1,1)*dgmetds(2,3))
       cm(4,4,1,3)=9*gmet(1,2)**2*gmet(3,3)*dgmetds(1,1)+gmet(1,3)**2*(9*gmet(2,2)&
&       *dgmetds(1,1)+6*gmet(1,2)*dgmetds(1,2))+gmet(1,1)*(9*gmet(2,3)&
&       **2*dgmetds(1,1)+15*gmet(2,2)*gmet(3,3)*dgmetds(1,1)+24*gmet(1,2)&
&       *gmet(3,3)*dgmetds(1,2)+6*gmet(1,2)*gmet(2,3)*dgmetds(1,3))-6*gmet(1,1)&
&       **2*gmet(2,3)*dgmetds(2,3)+gmet(1,3)*(6*gmet(1,2)**2*dgmetds(1,3)&
&       +gmet(1,1)*(6*gmet(2,3)*dgmetds(1,2)+24*gmet(2,2)*dgmetds(1,3))&
&       +gmet(1,2)*(-6*gmet(2,3)*dgmetds(1,1)+18*gmet(1,1)*dgmetds(2,3)))
       cm(5,4,1,3)=gmet(1,2)*(3*gmet(1,3)**2*dgmetds(1,1)+12*gmet(1,1)&
&       *gmet(3,3)*dgmetds(1,1)+6*gmet(1,1)*gmet(1,3)*dgmetds(1,3))+gmet(1,1)&
&       *(6*gmet(1,3)**2*dgmetds(1,2)+12*gmet(1,1)*(gmet(3,3)*dgmetds(1,2)&
&       +gmet(2,3)*dgmetds(1,3))+gmet(1,3)*(3*gmet(2,3)*dgmetds(1,1)&
&       +6*gmet(1,1)*dgmetds(2,3)))
       cm(6,4,1,3)=gmet(1,2)**2*(3*gmet(1,3)*dgmetds(1,1)+6*gmet(1,1)&
&       *dgmetds(1,3))+12*gmet(1,1)*(gmet(1,3)*gmet(2,2)*dgmetds(1,1)&
&       +gmet(1,1)*(gmet(2,3)*dgmetds(1,2)+gmet(2,2)*dgmetds(1,3)))+gmet(1,1)&
&       *gmet(1,2)*(3*gmet(2,3)*dgmetds(1,1)+6*(gmet(1,3)*dgmetds(1,2)&
&       +gmet(1,1)*dgmetds(2,3)))
       cm(7,4,1,3)=gmet(1,3)*gmet(2,2)*(-3*gmet(2,2)*dgmetds(1,1)-6*gmet(1,2)&
&       *dgmetds(1,2))-3*gmet(1,1)*gmet(2,2)*(gmet(2,3)*dgmetds(1,2)&
&       +gmet(2,2)*dgmetds(1,3))+gmet(1,2)**2*(15*gmet(2,3)*dgmetds(1,2)&
&       +9*gmet(2,2)*dgmetds(1,3))+5*gmet(1,2)**3*dgmetds(2,3)+gmet(1,2)&
&       *gmet(2,2)*(9*gmet(2,3)*dgmetds(1,1)-3*gmet(1,1)*dgmetds(2,3))
       cm(8,4,1,3)=-6*gmet(1,2)**2*gmet(3,3)*dgmetds(1,3)+gmet(1,1)*(-9*gmet(2,3)&
&       *gmet(3,3)*dgmetds(1,2)-6*gmet(2,3)**2*dgmetds(1,3)-3*gmet(2,2)&
&       *gmet(3,3)*dgmetds(1,3))+gmet(1,2)*gmet(3,3)*(3*gmet(2,3)*dgmetds(1,1)&
&       -3*gmet(1,1)*dgmetds(2,3))+gmet(1,3)*(3*gmet(2,3)**2*dgmetds(1,1)&
&       +gmet(3,3)*(12*gmet(2,2)*dgmetds(1,1)+24*gmet(1,2)*dgmetds(1,2))&
&       +gmet(2,3)*(18*gmet(1,2)*dgmetds(1,3)-6*gmet(1,1)*dgmetds(2,3)))&
&       +gmet(1,3)**2*(3*gmet(2,3)*dgmetds(1,2)+15*(gmet(2,2)*dgmetds(1,3)&
&       +gmet(1,2)*dgmetds(2,3)))
       cm(9,4,1,3)=-6*gmet(1,3)**2*gmet(2,2)*dgmetds(1,2)+gmet(1,2)**2*(15*gmet(3,3)&
&       *dgmetds(1,2)+3*gmet(2,3)*dgmetds(1,3))+gmet(1,1)*(-6*gmet(2,3)&
&       **2*dgmetds(1,2)-3*gmet(2,2)*gmet(3,3)*dgmetds(1,2)-9*gmet(2,2)&
&       *gmet(2,3)*dgmetds(1,3))+gmet(1,2)*(3*gmet(2,3)**2*dgmetds(1,1)&
&       +12*gmet(2,2)*gmet(3,3)*dgmetds(1,1)-6*gmet(1,1)*gmet(2,3)*dgmetds(2,3))&
&       +gmet(1,3)*(gmet(2,2)*(3*gmet(2,3)*dgmetds(1,1)+24*gmet(1,2)&
&       *dgmetds(1,3)-3*gmet(1,1)*dgmetds(2,3))+gmet(1,2)*(18*gmet(2,3)&
&       *dgmetds(1,2)+15*gmet(1,2)*dgmetds(2,3)))
       cm(10,4,1,3)=gmet(1,3)**2*(9*gmet(3,3)*dgmetds(1,2)+15*gmet(2,3)&
&       *dgmetds(1,3))-3*gmet(3,3)*(gmet(1,2)*gmet(3,3)*dgmetds(1,1)&
&       +gmet(1,1)*(gmet(3,3)*dgmetds(1,2)+gmet(2,3)*dgmetds(1,3)))+5*gmet(1,3)&
&       **3*dgmetds(2,3)+gmet(1,3)*gmet(3,3)*(9*gmet(2,3)*dgmetds(1,1)&
&       -6*gmet(1,2)*dgmetds(1,3)-3*gmet(1,1)*dgmetds(2,3))
       cm(1,5,1,3)=gmet(1,1)**2*(3*gmet(1,3)*dgmetds(1,1)+2*gmet(1,1)&
&       *dgmetds(1,3))
       cm(2,5,1,3)=12*gmet(1,1)*gmet(1,2)*gmet(2,3)*dgmetds(1,1)+gmet(1,1)&
&       *gmet(2,2)*(-4.5d0*gmet(1,3)*dgmetds(1,1)-3*gmet(1,1)*dgmetds(1,3))&
&       +gmet(1,2)**2*(1.5d0*gmet(1,3)*dgmetds(1,1)+9*gmet(1,1)*dgmetds(1,3))
       cm(3,5,1,3)=1.5d0*gmet(1,3)**3*dgmetds(1,1)+7.5d0*gmet(1,1)*gmet(1,3)&
&       *gmet(3,3)*dgmetds(1,1)+9*gmet(1,1)*gmet(1,3)**2*dgmetds(1,3)&
&       -3*gmet(1,1)**2*gmet(3,3)*dgmetds(1,3)
       cm(4,5,1,3)=gmet(1,1)*gmet(2,3)*(3*gmet(1,3)*dgmetds(1,1)-6*gmet(1,1)&
&       *dgmetds(1,3))+gmet(1,2)*(3*gmet(1,3)**2*dgmetds(1,1)+12*gmet(1,1)&
&       *gmet(3,3)*dgmetds(1,1)+18*gmet(1,1)*gmet(1,3)*dgmetds(1,3))
       cm(5,5,1,3)=gmet(1,1)*(3*gmet(1,3)**2*dgmetds(1,1)+6*gmet(1,1)&
&       *gmet(3,3)*dgmetds(1,1)+6*gmet(1,1)*gmet(1,3)*dgmetds(1,3))
       cm(6,5,1,3)=gmet(1,1)*(6*gmet(1,1)*gmet(2,3)*dgmetds(1,1)+gmet(1,2)&
&       *(3*gmet(1,3)*dgmetds(1,1)+6*gmet(1,1)*dgmetds(1,3)))
       cm(7,5,1,3)=7.5d0*gmet(1,2)**2*gmet(2,3)*dgmetds(1,1)-1.5d0*gmet(1,1)&
&       *gmet(2,2)*gmet(2,3)*dgmetds(1,1)+5*gmet(1,2)**3*dgmetds(1,3)&
&       -3*gmet(1,2)*gmet(2,2)*(gmet(1,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,3))
       cm(8,5,1,3)=gmet(1,1)*gmet(3,3)*(-4.5d0*gmet(2,3)*dgmetds(1,1)&
&       -3*gmet(1,2)*dgmetds(1,3))+gmet(1,3)**2*(1.5d0*gmet(2,3)*dgmetds(1,1)&
&       +15*gmet(1,2)*dgmetds(1,3))+gmet(1,3)*(12*gmet(1,2)*gmet(3,3)&
&       *dgmetds(1,1)-6*gmet(1,1)*gmet(2,3)*dgmetds(1,3))
       cm(9,5,1,3)=(12*(-12*gmet(1,3)**2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)+30*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(-12*gmet(2,3)&
&       **2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+8*(90*gmet(1,2)**2*gmet(1,3)&
&       -18*gmet(1,1)*gmet(1,3)*gmet(2,2)-36*gmet(1,1)*gmet(1,2)*gmet(2,3))&
&       *dgmetds(1,3))/48.d0
       cm(10,5,1,3)=4.5d0*gmet(1,3)**2*gmet(3,3)*dgmetds(1,1)-1.5d0*gmet(1,1)&
&       *gmet(3,3)**2*dgmetds(1,1)+5*gmet(1,3)**3*dgmetds(1,3)-3*gmet(1,1)&
&       *gmet(1,3)*gmet(3,3)*dgmetds(1,3)
       cm(1,6,1,3)=gmet(1,1)**2*(3*gmet(1,2)*dgmetds(1,1)+2*gmet(1,1)&
&       *dgmetds(1,2))
       cm(2,6,1,3)=1.5d0*gmet(1,2)**3*dgmetds(1,1)+7.5d0*gmet(1,1)*gmet(1,2)&
&       *gmet(2,2)*dgmetds(1,1)+9*gmet(1,1)*gmet(1,2)**2*dgmetds(1,2)&
&       -3*gmet(1,1)**2*gmet(2,2)*dgmetds(1,2)
       cm(3,6,1,3)=gmet(1,2)*(1.5d0*gmet(1,3)**2-4.5d0*gmet(1,1)*gmet(3,3))&
&       *dgmetds(1,1)+gmet(1,1)*(12*gmet(1,3)*gmet(2,3)*dgmetds(1,1)&
&       +9*gmet(1,3)**2*dgmetds(1,2)-3*gmet(1,1)*gmet(3,3)*dgmetds(1,2))
       cm(4,6,1,3)=3*gmet(1,2)**2*gmet(1,3)*dgmetds(1,1)+gmet(1,1)*gmet(1,2)&
&       *(3*gmet(2,3)*dgmetds(1,1)+18*gmet(1,3)*dgmetds(1,2))+gmet(1,1)&
&       *(12*gmet(1,3)*gmet(2,2)*dgmetds(1,1)-6*gmet(1,1)*gmet(2,3)*dgmetds(1,2))
       cm(5,6,1,3)=gmet(1,1)*(3*gmet(1,2)*gmet(1,3)*dgmetds(1,1)+6*gmet(1,1)&
&       *(gmet(2,3)*dgmetds(1,1)+gmet(1,3)*dgmetds(1,2)))
       cm(6,6,1,3)=gmet(1,1)*(3*gmet(1,2)**2*dgmetds(1,1)+6*gmet(1,1)&
&       *gmet(2,2)*dgmetds(1,1)+6*gmet(1,1)*gmet(1,2)*dgmetds(1,2))
       cm(7,6,1,3)=4.5d0*gmet(1,2)**2*gmet(2,2)*dgmetds(1,1)-1.5d0*gmet(1,1)&
&       *gmet(2,2)**2*dgmetds(1,1)+5*gmet(1,2)**3*dgmetds(1,2)-3*gmet(1,1)&
&       *gmet(1,2)*gmet(2,2)*dgmetds(1,2)
       cm(8,6,1,3)=(12*(30*gmet(1,3)**2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(-12*gmet(2,3)&
&       **2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+8*(-36*gmet(1,1)*gmet(1,3)&
&       *gmet(2,3)+gmet(1,2)*(90*gmet(1,3)**2-18*gmet(1,1)*gmet(3,3)))&
&       *dgmetds(1,2))/48.d0
       cm(9,6,1,3)=gmet(1,1)*gmet(2,2)*(-4.5d0*gmet(2,3)*dgmetds(1,1)&
&       -3*gmet(1,3)*dgmetds(1,2))+gmet(1,2)**2*(1.5d0*gmet(2,3)*dgmetds(1,1)&
&       +15*gmet(1,3)*dgmetds(1,2))+gmet(1,2)*(12*gmet(1,3)*gmet(2,2)&
&       *dgmetds(1,1)-6*gmet(1,1)*gmet(2,3)*dgmetds(1,2))
       cm(10,6,1,3)=7.5d0*gmet(1,3)**2*gmet(2,3)*dgmetds(1,1)-1.5d0*gmet(1,1)&
&       *gmet(2,3)*gmet(3,3)*dgmetds(1,1)+5*gmet(1,3)**3*dgmetds(1,2)&
&       -3*gmet(1,3)*gmet(3,3)*(gmet(1,2)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,2))
       cm(1,7,1,3)=2.5d0*gmet(1,2)**3*dgmetds(1,1)+9*gmet(1,1)*gmet(1,2)&
&       **2*dgmetds(1,2)-3*gmet(1,1)**2*gmet(2,2)*dgmetds(1,2)+gmet(1,1)&
&       *gmet(1,2)*(-1.5d0*gmet(2,2)*dgmetds(1,1)+3*gmet(1,1)*dgmetds(2,2))
       cm(2,7,1,3)=6*gmet(1,2)**2*gmet(2,2)*dgmetds(1,2)+12*gmet(1,1)&
&       *gmet(2,2)**2*dgmetds(1,2)+1.5d0*gmet(1,2)**3*dgmetds(2,2)+gmet(1,2)&
&       *gmet(2,2)*(3*gmet(2,2)*dgmetds(1,1)+7.5d0*gmet(1,1)*dgmetds(2,2))
       cm(3,7,1,3)=-6*gmet(1,2)**2*gmet(3,3)*dgmetds(1,2)+gmet(1,1)*(15*gmet(2,3)&
&       **2-3*gmet(2,2)*gmet(3,3))*dgmetds(1,2)+gmet(1,3)*gmet(2,3)*(-3*gmet(2,2)&
&       *dgmetds(1,1)+18*gmet(1,2)*dgmetds(1,2)+12*gmet(1,1)*dgmetds(2,2))&
&       +gmet(1,3)**2*(-6*gmet(2,2)*dgmetds(1,2)+1.5d0*gmet(1,2)*dgmetds(2,2))&
&       +gmet(1,2)*(7.5d0*gmet(2,3)**2*dgmetds(1,1)-1.5d0*gmet(2,2)*gmet(3,3)&
&       *dgmetds(1,1)-4.5d0*gmet(1,1)*gmet(3,3)*dgmetds(2,2))
       cm(4,7,1,3)=gmet(2,3)*(6*gmet(1,2)**2*dgmetds(1,2)+24*gmet(1,1)&
&       *gmet(2,2)*dgmetds(1,2)+gmet(1,2)*(9*gmet(2,2)*dgmetds(1,1)+3*gmet(1,1)&
&       *dgmetds(2,2)))+gmet(1,3)*(-3*gmet(2,2)**2*dgmetds(1,1)+3*gmet(1,2)&
&       **2*dgmetds(2,2)+gmet(2,2)*(6*gmet(1,2)*dgmetds(1,2)+12*gmet(1,1)&
&       *dgmetds(2,2)))
       cm(5,7,1,3)=gmet(1,2)**2*(7.5d0*gmet(2,3)*dgmetds(1,1)+3*gmet(1,3)&
&       *dgmetds(1,2))+gmet(1,1)*(gmet(2,2)*(-1.5d0*gmet(2,3)*dgmetds(1,1)&
&       -9*gmet(1,3)*dgmetds(1,2))+6*gmet(1,1)*gmet(2,3)*dgmetds(2,2))&
&       +gmet(1,2)*(24*gmet(1,1)*gmet(2,3)*dgmetds(1,2)+gmet(1,3)*(-3*gmet(2,2)&
&       *dgmetds(1,1)+3*gmet(1,1)*dgmetds(2,2)))
       cm(6,7,1,3)=3*gmet(1,2)**3*dgmetds(1,2)+15*gmet(1,1)*gmet(1,2)&
&       *gmet(2,2)*dgmetds(1,2)+gmet(1,2)**2*(4.5d0*gmet(2,2)*dgmetds(1,1)&
&       +3*gmet(1,1)*dgmetds(2,2))+gmet(1,1)*gmet(2,2)*(-1.5d0*gmet(2,2)&
&       *dgmetds(1,1)+6*gmet(1,1)*dgmetds(2,2))
       cm(7,7,1,3)=gmet(2,2)*(1*gmet(2,2)**2*dgmetds(1,1)+4.5d0*gmet(1,2)&
&       **2*dgmetds(2,2)+gmet(2,2)*(6*gmet(1,2)*dgmetds(1,2)-1.5d0*gmet(1,1)&
&       *dgmetds(2,2)))
       cm(8,7,1,3)=(2*gmet(2,2)*(54*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3))&
&       *dgmetds(1,1)+12*(48*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)&
&       *(6*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+6*(30*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(-12*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(2,2))&
&       /24.d0
       cm(9,7,1,3)=gmet(2,2)**2*(3*gmet(2,3)*dgmetds(1,1)+12*gmet(1,3)&
&       *dgmetds(1,2))+1.5d0*gmet(1,2)**2*gmet(2,3)*dgmetds(2,2)+gmet(2,2)&
&       *(-4.5d0*gmet(1,1)*gmet(2,3)*dgmetds(2,2)+gmet(1,2)*(6*gmet(2,3)&
&       *dgmetds(1,2)+12*gmet(1,3)*dgmetds(2,2)))
       cm(10,7,1,3)=(2*(90*gmet(2,3)**3-54*gmet(2,2)*gmet(2,3)*gmet(3,3))&
&       *dgmetds(1,1)+12*(-36*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)&
&       *(90*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+6*(90*gmet(1,3)&
&       **2*gmet(2,3)-36*gmet(1,2)*gmet(1,3)*gmet(3,3)-18*gmet(1,1)*gmet(2,3)&
&       *gmet(3,3))*dgmetds(2,2))/72.d0
       cm(1,8,1,3)=(gmet(1,1)*(216*gmet(1,3)**2*dgmetds(1,2)+gmet(1,1)&
&       *(-72*gmet(3,3)*dgmetds(1,2)-144*gmet(2,3)*dgmetds(1,3))+gmet(1,3)&
&       *(-72*gmet(2,3)*dgmetds(1,1)+144*gmet(1,1)*dgmetds(2,3)))+gmet(1,2)&
&       *(180*gmet(1,3)**2*dgmetds(1,1)+432*gmet(1,1)*gmet(1,3)*dgmetds(1,3)&
&       +gmet(1,1)*(-36*gmet(3,3)*dgmetds(1,1)+72*gmet(1,1)*dgmetds(3,3))))&
&       /24.d0
       cm(2,8,1,3)=-6*gmet(1,3)**2*gmet(2,2)*dgmetds(1,2)+gmet(1,2)**2*(-6*gmet(3,3)&
&       *dgmetds(1,2)+6*gmet(2,3)*dgmetds(1,3))+gmet(1,1)*(15*gmet(2,3)&
&       **2*dgmetds(1,2)-3*gmet(2,2)*gmet(3,3)*dgmetds(1,2)+24*gmet(2,2)&
&       *gmet(2,3)*dgmetds(1,3))+gmet(1,3)*(gmet(2,2)*(12*gmet(2,3)*dgmetds(1,1)&
&       +6*gmet(1,2)*dgmetds(1,3)-9*gmet(1,1)*dgmetds(2,3))+gmet(1,2)&
&       *(18*gmet(2,3)*dgmetds(1,2)+3*gmet(1,2)*dgmetds(2,3)))+1.5d0*gmet(1,2)&
&       **3*dgmetds(3,3)+gmet(1,2)*(1.5d0*gmet(2,3)**2*dgmetds(1,1)-4.5d0*gmet(2,2)&
&       *gmet(3,3)*dgmetds(1,1)+24*gmet(1,1)*gmet(2,3)*dgmetds(2,3)+7.5d0*gmet(1,1)&
&       *gmet(2,2)*dgmetds(3,3))
       cm(3,8,1,3)=3*gmet(1,3)**3*dgmetds(2,3)+gmet(1,3)**2*(6*gmet(3,3)&
&       *dgmetds(1,2)+6*gmet(2,3)*dgmetds(1,3)+1.5d0*gmet(1,2)*dgmetds(3,3))&
&       +gmet(3,3)*(gmet(1,1)*(12*gmet(3,3)*dgmetds(1,2)+24*gmet(2,3)&
&       *dgmetds(1,3))+gmet(1,2)*(6*gmet(3,3)*dgmetds(1,1)-4.5d0*gmet(1,1)&
&       *dgmetds(3,3)))+gmet(1,3)*(gmet(3,3)*(6*gmet(1,2)*dgmetds(1,3)&
&       +15*gmet(1,1)*dgmetds(2,3))+gmet(2,3)*(3*gmet(3,3)*dgmetds(1,1)&
&       +12*gmet(1,1)*dgmetds(3,3)))
       cm(4,8,1,3)=18*gmet(1,2)**2*gmet(3,3)*dgmetds(1,3)+gmet(1,1)*(24*gmet(2,3)&
&       *gmet(3,3)*dgmetds(1,2)+18*gmet(2,3)**2*dgmetds(1,3)+30*gmet(2,2)&
&       *gmet(3,3)*dgmetds(1,3))+gmet(1,3)**2*(6*gmet(2,3)*dgmetds(1,2)&
&       +18*gmet(2,2)*dgmetds(1,3)+6*gmet(1,2)*dgmetds(2,3))+gmet(1,2)&
&       *(24*gmet(1,1)*gmet(3,3)*dgmetds(2,3)+3*gmet(2,3)*(gmet(3,3)&
&       *dgmetds(1,1)+gmet(1,1)*dgmetds(3,3)))+gmet(1,3)*(3*gmet(2,3)&
&       **2*dgmetds(1,1)+gmet(2,3)*(-12*gmet(1,2)*dgmetds(1,3)+6*gmet(1,1)&
&       *dgmetds(2,3))+12*gmet(2,2)*(gmet(3,3)*dgmetds(1,1)+gmet(1,1)&
&       *dgmetds(3,3))+gmet(1,2)*(6*gmet(3,3)*dgmetds(1,2)+3*gmet(1,2)&
&       *dgmetds(3,3)))
       cm(5,8,1,3)=3*gmet(1,3)**3*dgmetds(1,2)+gmet(1,3)**2*(1.5d0*gmet(2,3)&
&       *dgmetds(1,1)+6*(gmet(1,2)*dgmetds(1,3)+gmet(1,1)*dgmetds(2,3)))&
&       +gmet(1,3)*(gmet(1,1)*(15*gmet(3,3)*dgmetds(1,2)+6*gmet(2,3)&
&       *dgmetds(1,3))+gmet(1,2)*(12*gmet(3,3)*dgmetds(1,1)+3*gmet(1,1)&
&       *dgmetds(3,3)))+gmet(1,1)*(gmet(3,3)*(24*gmet(1,2)*dgmetds(1,3)&
&       +12*gmet(1,1)*dgmetds(2,3))+gmet(2,3)*(-4.5d0*gmet(3,3)*dgmetds(1,1)&
&       +6*gmet(1,1)*dgmetds(3,3)))
       cm(6,8,1,3)=-3*gmet(1,2)**2*gmet(3,3)*dgmetds(1,1)+gmet(1,3)**2*(7.5d0*gmet(2,2)&
&       *dgmetds(1,1)+3*gmet(1,2)*dgmetds(1,2))+gmet(1,3)*(6*gmet(1,2)&
&       **2*dgmetds(1,3)+24*gmet(1,1)*(gmet(2,3)*dgmetds(1,2)+gmet(2,2)&
&       *dgmetds(1,3))+gmet(1,2)*(9*gmet(2,3)*dgmetds(1,1)+6*gmet(1,1)&
&       *dgmetds(2,3)))+gmet(1,1)*(-3*gmet(2,3)**2*dgmetds(1,1)-1.5d0*gmet(2,2)&
&       *gmet(3,3)*dgmetds(1,1)-9*gmet(1,2)*gmet(3,3)*dgmetds(1,2)+6*gmet(1,2)&
&       *gmet(2,3)*dgmetds(1,3)+3*gmet(1,2)**2*dgmetds(3,3))+gmet(1,1)&
&       **2*(12*gmet(2,3)*dgmetds(2,3)+6*gmet(2,2)*dgmetds(3,3))
       cm(7,8,1,3)=15*gmet(1,2)*gmet(2,3)*(gmet(2,3)*dgmetds(1,2)+gmet(1,2)&
&       *dgmetds(2,3))+gmet(2,2)**2*(-1.5d0*gmet(3,3)*dgmetds(1,1)-6*gmet(1,3)&
&       *dgmetds(1,3)-1.5d0*gmet(1,1)*dgmetds(3,3))+gmet(2,2)*(4.5d0*gmet(2,3)&
&       **2*dgmetds(1,1)+gmet(2,3)*(-6*gmet(1,3)*dgmetds(1,2)+18*gmet(1,2)&
&       *dgmetds(1,3)-3*gmet(1,1)*dgmetds(2,3))+gmet(1,2)*(-3*gmet(3,3)&
&       *dgmetds(1,2)-6*gmet(1,3)*dgmetds(2,3)+4.5d0*gmet(1,2)*dgmetds(3,3)))
       cm(8,8,1,3)=gmet(2,3)**2*(3*gmet(3,3)*dgmetds(1,1)+6*gmet(1,3)&
&       *dgmetds(1,3)-3*gmet(1,1)*dgmetds(3,3))+gmet(1,2)*gmet(3,3)*(12*gmet(3,3)&
&       *dgmetds(1,2)+24*gmet(1,3)*dgmetds(2,3)-3*gmet(1,2)*dgmetds(3,3))&
&       +gmet(2,2)*(6*gmet(3,3)**2*dgmetds(1,1)+24*gmet(1,3)*gmet(3,3)&
&       *dgmetds(1,3)+7.5d0*gmet(1,3)**2*dgmetds(3,3)-1.5d0*gmet(1,1)&
&       *gmet(3,3)*dgmetds(3,3))+gmet(2,3)*(3*gmet(1,3)**2*dgmetds(2,3)&
&       +gmet(3,3)*(6*gmet(1,2)*dgmetds(1,3)-9*gmet(1,1)*dgmetds(2,3))&
&       +gmet(1,3)*(6*gmet(3,3)*dgmetds(1,2)+9*gmet(1,2)*dgmetds(3,3)))
       cm(9,8,1,3)=1.5d0*gmet(2,3)**3*dgmetds(1,1)-6*gmet(1,3)**2*gmet(2,2)&
&       *dgmetds(2,3)+gmet(2,3)**2*(3*gmet(1,3)*dgmetds(1,2)+6*gmet(1,2)&
&       *dgmetds(1,3)-6*gmet(1,1)*dgmetds(2,3))+gmet(3,3)*(24*gmet(1,2)&
&       *gmet(2,2)*dgmetds(1,3)+15*gmet(1,2)**2*dgmetds(2,3)-3*gmet(1,1)&
&       *gmet(2,2)*dgmetds(2,3))+gmet(1,3)*gmet(2,2)*(-9*gmet(3,3)*dgmetds(1,2)&
&       +12*gmet(1,2)*dgmetds(3,3))+gmet(2,3)*(gmet(2,2)*(7.5d0*gmet(3,3)&
&       *dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,3)-4.5d0*gmet(1,1)*dgmetds(3,3))&
&       +gmet(1,2)*(24*gmet(3,3)*dgmetds(1,2)+18*gmet(1,3)*dgmetds(2,3)&
&       +1.5d0*gmet(1,2)*dgmetds(3,3)))
       cm(10,8,1,3)=gmet(2,3)*(3*gmet(3,3)**2*dgmetds(1,1)+7.5d0*gmet(1,3)&
&       **2*dgmetds(3,3)+gmet(3,3)*(18*gmet(1,3)*dgmetds(1,3)-1.5d0*gmet(1,1)&
&       *dgmetds(3,3)))+gmet(3,3)*(9*gmet(1,3)**2*dgmetds(2,3)+gmet(3,3)&
&       *(-6*gmet(1,2)*dgmetds(1,3)-3*gmet(1,1)*dgmetds(2,3))+gmet(1,3)&
&       *(6*gmet(3,3)*dgmetds(1,2)-3*gmet(1,2)*dgmetds(3,3)))
       cm(1,9,1,3)=(gmet(1,2)**2*(180*gmet(1,3)*dgmetds(1,1)+216*gmet(1,1)&
&       *dgmetds(1,3))+gmet(1,1)*(gmet(1,1)*(-144*gmet(2,3)*dgmetds(1,2)&
&       -72*gmet(2,2)*dgmetds(1,3))+gmet(1,3)*(-36*gmet(2,2)*dgmetds(1,1)&
&       +72*gmet(1,1)*dgmetds(2,2)))+gmet(1,1)*gmet(1,2)*(-72*gmet(2,3)&
&       *dgmetds(1,1)+432*gmet(1,3)*dgmetds(1,2)+144*gmet(1,1)*dgmetds(2,3)))&
&       /24.d0
       cm(2,9,1,3)=6*gmet(1,2)**2*(gmet(2,3)*dgmetds(1,2)+gmet(2,2)*dgmetds(1,3))&
&       +gmet(1,1)*gmet(2,2)*(24*gmet(2,3)*dgmetds(1,2)+12*gmet(2,2)&
&       *dgmetds(1,3))+gmet(1,3)*(6*gmet(2,2)**2*dgmetds(1,1)+1.5d0*gmet(1,2)&
&       **2*dgmetds(2,2)+gmet(2,2)*(6*gmet(1,2)*dgmetds(1,2)-4.5d0*gmet(1,1)&
&       *dgmetds(2,2)))+3*gmet(1,2)**3*dgmetds(2,3)+gmet(1,2)*(3*gmet(2,2)&
&       *gmet(2,3)*dgmetds(1,1)+12*gmet(1,1)*gmet(2,3)*dgmetds(2,2)+15*gmet(1,1)&
&       *gmet(2,2)*dgmetds(2,3))
       cm(3,9,1,3)=-6*gmet(1,2)**2*gmet(3,3)*dgmetds(1,3)+gmet(1,1)*(24*gmet(2,3)&
&       *gmet(3,3)*dgmetds(1,2)+15*gmet(2,3)**2*dgmetds(1,3)-3*gmet(2,2)&
&       *gmet(3,3)*dgmetds(1,3))+1.5d0*gmet(1,3)**3*dgmetds(2,2)+gmet(1,2)&
&       *gmet(3,3)*(12*gmet(2,3)*dgmetds(1,1)-9*gmet(1,1)*dgmetds(2,3))&
&       +gmet(1,3)**2*(6*gmet(2,3)*dgmetds(1,2)-6*gmet(2,2)*dgmetds(1,3)&
&       +3*gmet(1,2)*dgmetds(2,3))+gmet(1,3)*(1.5d0*gmet(2,3)**2*dgmetds(1,1)&
&       +gmet(3,3)*(-4.5d0*gmet(2,2)*dgmetds(1,1)+6*gmet(1,2)*dgmetds(1,2)&
&       +7.5d0*gmet(1,1)*dgmetds(2,2))+gmet(2,3)*(18*gmet(1,2)*dgmetds(1,3)&
&       +24*gmet(1,1)*dgmetds(2,3)))
       cm(4,9,1,3)=gmet(1,2)**2*(18*gmet(3,3)*dgmetds(1,2)+6*gmet(2,3)&
&       *dgmetds(1,3))+gmet(1,1)*(18*gmet(2,3)**2*dgmetds(1,2)+30*gmet(2,2)&
&       *gmet(3,3)*dgmetds(1,2)+24*gmet(2,2)*gmet(2,3)*dgmetds(1,3))&
&       +gmet(1,3)**2*(18*gmet(2,2)*dgmetds(1,2)+3*gmet(1,2)*dgmetds(2,2))&
&       +gmet(1,2)*(3*gmet(2,3)**2*dgmetds(1,1)+12*gmet(3,3)*(gmet(2,2)&
&       *dgmetds(1,1)+gmet(1,1)*dgmetds(2,2))+6*gmet(1,1)*gmet(2,3)*dgmetds(2,3))&
&       +gmet(1,3)*(-12*gmet(1,2)*gmet(2,3)*dgmetds(1,2)+3*gmet(1,1)&
&       *gmet(2,3)*dgmetds(2,2)+6*gmet(1,2)**2*dgmetds(2,3)+gmet(2,2)&
&       *(3*gmet(2,3)*dgmetds(1,1)+6*gmet(1,2)*dgmetds(1,3)+24*gmet(1,1)&
&       *dgmetds(2,3)))
       cm(5,9,1,3)=7.5d0*gmet(1,2)**2*gmet(3,3)*dgmetds(1,1)+gmet(1,1)&
&       *(-3*gmet(2,3)**2*dgmetds(1,1)-1.5d0*gmet(2,2)*gmet(3,3)*dgmetds(1,1)&
&       +24*gmet(1,2)*gmet(3,3)*dgmetds(1,2)+24*gmet(1,2)*gmet(2,3)*dgmetds(1,3))&
&       +gmet(1,3)**2*(-3*gmet(2,2)*dgmetds(1,1)+6*gmet(1,2)*dgmetds(1,2)&
&       +3*gmet(1,1)*dgmetds(2,2))+gmet(1,1)**2*(6*gmet(3,3)*dgmetds(2,2)&
&       +12*gmet(2,3)*dgmetds(2,3))+gmet(1,3)*(3*gmet(1,2)**2*dgmetds(1,3)&
&       +gmet(1,1)*(6*gmet(2,3)*dgmetds(1,2)-9*gmet(2,2)*dgmetds(1,3))&
&       +gmet(1,2)*(9*gmet(2,3)*dgmetds(1,1)+6*gmet(1,1)*dgmetds(2,3)))
       cm(6,9,1,3)=3*gmet(1,2)**3*dgmetds(1,3)+gmet(1,2)*(gmet(1,1)*(6*gmet(2,3)&
&       *dgmetds(1,2)+15*gmet(2,2)*dgmetds(1,3))+gmet(1,3)*(12*gmet(2,2)&
&       *dgmetds(1,1)+3*gmet(1,1)*dgmetds(2,2)))+gmet(1,2)**2*(1.5d0*gmet(2,3)&
&       *dgmetds(1,1)+6*(gmet(1,3)*dgmetds(1,2)+gmet(1,1)*dgmetds(2,3)))&
&       +gmet(1,1)*(6*gmet(1,1)*gmet(2,3)*dgmetds(2,2)+gmet(2,2)*(-4.5d0*gmet(2,3)&
&       *dgmetds(1,1)+24*gmet(1,3)*dgmetds(1,2)+12*gmet(1,1)*dgmetds(2,3)))
       cm(7,9,1,3)=7.5d0*gmet(1,2)**2*gmet(2,3)*dgmetds(2,2)+gmet(2,2)&
&       **2*(3*gmet(2,3)*dgmetds(1,1)-6*gmet(1,3)*dgmetds(1,2)+6*gmet(1,2)&
&       *dgmetds(1,3)-3*gmet(1,1)*dgmetds(2,3))+gmet(2,2)*(-1.5d0*gmet(1,1)&
&       *gmet(2,3)*dgmetds(2,2)+gmet(1,2)*(18*gmet(2,3)*dgmetds(1,2)&
&       -3*gmet(1,3)*dgmetds(2,2))+9*gmet(1,2)**2*dgmetds(2,3))
       cm(8,9,1,3)=1.5d0*gmet(2,3)**3*dgmetds(1,1)+gmet(1,3)*gmet(3,3)&
&       *(24*gmet(2,2)*dgmetds(1,2)+12*gmet(1,2)*dgmetds(2,2))+15*gmet(1,3)&
&       **2*gmet(2,2)*dgmetds(2,3)+gmet(2,3)**2*(6*gmet(1,3)*dgmetds(1,2)&
&       +3*gmet(1,2)*dgmetds(1,3)-6*gmet(1,1)*dgmetds(2,3))+gmet(3,3)&
&       *(-9*gmet(1,2)*gmet(2,2)*dgmetds(1,3)-6*gmet(1,2)**2*dgmetds(2,3)&
&       -3*gmet(1,1)*gmet(2,2)*dgmetds(2,3))+gmet(2,3)*(gmet(2,2)*(7.5d0*gmet(3,3)&
&       *dgmetds(1,1)+24*gmet(1,3)*dgmetds(1,3))+(1.5d0*gmet(1,3)**2-4.5d0*gmet(1,1)&
&       *gmet(3,3))*dgmetds(2,2)+gmet(1,2)*(6*gmet(3,3)*dgmetds(1,2)&
&       +18*gmet(1,3)*dgmetds(2,3)))
       cm(9,9,1,3)=gmet(2,2)**2*(6*gmet(3,3)*dgmetds(1,1)+12*gmet(1,3)&
&       *dgmetds(1,3))-3*gmet(1,1)*gmet(2,3)**2*dgmetds(2,2)+gmet(1,2)&
&       *gmet(2,3)*(6*gmet(2,3)*dgmetds(1,2)+9*gmet(1,3)*dgmetds(2,2))&
&       +gmet(1,2)**2*(7.5d0*gmet(3,3)*dgmetds(2,2)+3*gmet(2,3)*dgmetds(2,3))&
&       +gmet(2,2)*(3*gmet(2,3)**2*dgmetds(1,1)+(-3*gmet(1,3)**2-1.5d0*gmet(1,1)&
&       *gmet(3,3))*dgmetds(2,2)+gmet(2,3)*(6*gmet(1,3)*dgmetds(1,2)&
&       +6*gmet(1,2)*dgmetds(1,3)-9*gmet(1,1)*dgmetds(2,3))+24*gmet(1,2)&
&       *(gmet(3,3)*dgmetds(1,2)+gmet(1,3)*dgmetds(2,3)))
       cm(10,9,1,3)=gmet(2,3)**2*(4.5d0*gmet(3,3)*dgmetds(1,1)+15*gmet(1,3)&
&       *dgmetds(1,3))+gmet(2,3)*(18*gmet(1,3)*gmet(3,3)*dgmetds(1,2)&
&       +15*gmet(1,3)**2*dgmetds(2,3)+gmet(3,3)*(-6*gmet(1,2)*dgmetds(1,3)&
&       -3*gmet(1,1)*dgmetds(2,3)))+gmet(3,3)*(gmet(2,2)*(-1.5d0*gmet(3,3)&
&       *dgmetds(1,1)-3*gmet(1,3)*dgmetds(1,3))+(4.5d0*gmet(1,3)**2-1.5d0*gmet(1,1)&
&       *gmet(3,3))*dgmetds(2,2)-6*gmet(1,2)*(gmet(3,3)*dgmetds(1,2)&
&       +gmet(1,3)*dgmetds(2,3)))
       cm(1,10,1,3)=2.5d0*gmet(1,3)**3*dgmetds(1,1)+9*gmet(1,1)*gmet(1,3)&
&       **2*dgmetds(1,3)-3*gmet(1,1)**2*gmet(3,3)*dgmetds(1,3)+gmet(1,1)&
&       *gmet(1,3)*(-1.5d0*gmet(3,3)*dgmetds(1,1)+3*gmet(1,1)*dgmetds(3,3))
       cm(2,10,1,3)=-6*gmet(1,3)**2*gmet(2,2)*dgmetds(1,3)-6*gmet(1,2)&
&       **2*gmet(3,3)*dgmetds(1,3)+gmet(1,1)*(15*gmet(2,3)**2-3*gmet(2,2)&
&       *gmet(3,3))*dgmetds(1,3)+gmet(1,2)*gmet(2,3)*(-3*gmet(3,3)*dgmetds(1,1)&
&       +12*gmet(1,1)*dgmetds(3,3))+gmet(1,3)*(7.5d0*gmet(2,3)**2*dgmetds(1,1)&
&       -1.5d0*gmet(2,2)*gmet(3,3)*dgmetds(1,1)+18*gmet(1,2)*gmet(2,3)&
&       *dgmetds(1,3)+1.5d0*gmet(1,2)**2*dgmetds(3,3)-4.5d0*gmet(1,1)&
&       *gmet(2,2)*dgmetds(3,3))
       cm(3,10,1,3)=6*gmet(1,3)**2*gmet(3,3)*dgmetds(1,3)+12*gmet(1,1)&
&       *gmet(3,3)**2*dgmetds(1,3)+1.5d0*gmet(1,3)**3*dgmetds(3,3)+gmet(1,3)&
&       *gmet(3,3)*(3*gmet(3,3)*dgmetds(1,1)+7.5d0*gmet(1,1)*dgmetds(3,3))
       cm(4,10,1,3)=gmet(1,3)**2*(6*gmet(2,3)*dgmetds(1,3)+3*gmet(1,2)&
&       *dgmetds(3,3))+gmet(1,3)*(6*gmet(1,2)*gmet(3,3)*dgmetds(1,3)&
&       +gmet(2,3)*(9*gmet(3,3)*dgmetds(1,1)+3*gmet(1,1)*dgmetds(3,3)))&
&       +gmet(3,3)*(24*gmet(1,1)*gmet(2,3)*dgmetds(1,3)+gmet(1,2)*(-3*gmet(3,3)&
&       *dgmetds(1,1)+12*gmet(1,1)*dgmetds(3,3)))
       cm(5,10,1,3)=3*gmet(1,3)**3*dgmetds(1,3)+15*gmet(1,1)*gmet(1,3)&
&       *gmet(3,3)*dgmetds(1,3)+gmet(1,3)**2*(4.5d0*gmet(3,3)*dgmetds(1,1)&
&       +3*gmet(1,1)*dgmetds(3,3))+gmet(1,1)*gmet(3,3)*(-1.5d0*gmet(3,3)&
&       *dgmetds(1,1)+6*gmet(1,1)*dgmetds(3,3))
       cm(6,10,1,3)=gmet(1,3)**2*(7.5d0*gmet(2,3)*dgmetds(1,1)+3*gmet(1,2)&
&       *dgmetds(1,3))+gmet(1,3)*(24*gmet(1,1)*gmet(2,3)*dgmetds(1,3)&
&       +gmet(1,2)*(-3*gmet(3,3)*dgmetds(1,1)+3*gmet(1,1)*dgmetds(3,3)))&
&       +gmet(1,1)*(-9*gmet(1,2)*gmet(3,3)*dgmetds(1,3)+gmet(2,3)*(-1.5d0*gmet(3,3)&
&       *dgmetds(1,1)+6*gmet(1,1)*dgmetds(3,3)))
       cm(7,10,1,3)=(180*gmet(2,3)**3*dgmetds(1,1)+1080*gmet(1,2)*gmet(2,3)&
&       **2*dgmetds(1,3)-216*gmet(1,2)*gmet(2,2)*(gmet(3,3)*dgmetds(1,3)&
&       +gmet(1,3)*dgmetds(3,3))+gmet(2,3)*(540*gmet(1,2)**2*dgmetds(3,3)&
&       +gmet(2,2)*(-108*gmet(3,3)*dgmetds(1,1)-432*gmet(1,3)*dgmetds(1,3)&
&       -108*gmet(1,1)*dgmetds(3,3))))/72.d0
       cm(8,10,1,3)=12*gmet(1,2)*gmet(3,3)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)&
&       *dgmetds(3,3))+gmet(2,3)*(3*gmet(3,3)**2*dgmetds(1,1)+1.5d0*gmet(1,3)&
&       **2*dgmetds(3,3)+gmet(3,3)*(6*gmet(1,3)*dgmetds(1,3)-4.5d0*gmet(1,1)&
&       *dgmetds(3,3)))
       cm(9,10,1,3)=(2*gmet(3,3)*(54*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3))&
&       *dgmetds(1,1)+12*(48*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)&
&       *(6*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,3)+6*(-12*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)+30*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(-12*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(3,3))&
&       /24.d0
       cm(10,10,1,3)=gmet(3,3)*(1*gmet(3,3)**2*dgmetds(1,1)+4.5d0*gmet(1,3)&
&       **2*dgmetds(3,3)+gmet(3,3)*(6*gmet(1,3)*dgmetds(1,3)-1.5d0*gmet(1,1)&
&       *dgmetds(3,3)))
       cm(1,11,1,3)=5*gmet(1,2)**3*dgmetds(1,2)-3*gmet(1,1)*gmet(1,2)&
&       *gmet(2,2)*dgmetds(1,2)+4.5d0*gmet(1,1)*gmet(1,2)**2*dgmetds(2,2)&
&       -1.5d0*gmet(1,1)**2*gmet(2,2)*dgmetds(2,2)
       cm(2,11,1,3)=gmet(2,2)*(6*gmet(1,2)*gmet(2,2)*dgmetds(1,2)+3*gmet(1,2)&
&       **2*dgmetds(2,2)+6*gmet(1,1)*gmet(2,2)*dgmetds(2,2))
       cm(3,11,1,3)=(8*(-36*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+12*(-12*gmet(1,3)**2*gmet(2,2)&
&       +36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)&
&       *(30*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(2,2))/48.d0
       cm(4,11,1,3)=gmet(1,3)*gmet(2,2)*(-6*gmet(2,2)*dgmetds(1,2)+3*gmet(1,2)&
&       *dgmetds(2,2))+gmet(2,3)*(18*gmet(1,2)*gmet(2,2)*dgmetds(1,2)&
&       +3*gmet(1,2)**2*dgmetds(2,2)+12*gmet(1,1)*gmet(2,2)*dgmetds(2,2))
       cm(5,11,1,3)=(8*(-36*gmet(1,2)*gmet(1,3)*gmet(2,2)+90*gmet(1,2)&
&       **2*gmet(2,3)-18*gmet(1,1)*gmet(2,2)*gmet(2,3))*dgmetds(1,2)&
&       +12*(6*gmet(1,2)**2*gmet(1,3)-18*gmet(1,1)*gmet(1,3)*gmet(2,2)&
&       +48*gmet(1,1)*gmet(1,2)*gmet(2,3))*dgmetds(2,2))/48.d0
       cm(6,11,1,3)=9*gmet(1,2)**2*gmet(2,2)*dgmetds(1,2)-3*gmet(1,1)&
&       *gmet(2,2)**2*dgmetds(1,2)+1.5d0*gmet(1,2)**3*dgmetds(2,2)+7.5d0*gmet(1,1)&
&       *gmet(1,2)*gmet(2,2)*dgmetds(2,2)
       cm(7,11,1,3)=gmet(2,2)**2*(2*gmet(2,2)*dgmetds(1,2)+3*gmet(1,2)&
&       *dgmetds(2,2))
       cm(8,11,1,3)=-3*gmet(2,2)**2*gmet(3,3)*dgmetds(1,2)+1.5d0*gmet(1,2)&
&       *gmet(2,3)**2*dgmetds(2,2)+gmet(2,2)*(9*gmet(2,3)**2*dgmetds(1,2)&
&       +12*gmet(1,3)*gmet(2,3)*dgmetds(2,2)-4.5d0*gmet(1,2)*gmet(3,3)&
&       *dgmetds(2,2))
       cm(9,11,1,3)=gmet(2,2)*(3*gmet(1,2)*gmet(2,3)*dgmetds(2,2)+6*gmet(2,2)&
&       *(gmet(2,3)*dgmetds(1,2)+gmet(1,3)*dgmetds(2,2)))
       cm(10,11,1,3)=5*gmet(2,3)**3*dgmetds(1,2)+7.5d0*gmet(1,3)*gmet(2,3)&
&       **2*dgmetds(2,2)-1.5d0*gmet(1,3)*gmet(2,2)*gmet(3,3)*dgmetds(2,2)&
&       -3*gmet(2,3)*gmet(3,3)*(gmet(2,2)*dgmetds(1,2)+gmet(1,2)*dgmetds(2,2))
       cm(1,12,1,3)=(2*(-36*gmet(1,1)*gmet(1,3)*gmet(2,3)+gmet(1,2)*(90*gmet(1,3)&
&       **2-18*gmet(1,1)*gmet(3,3)))*dgmetds(1,2)+2*(90*gmet(1,2)**2*gmet(1,3)&
&       -18*gmet(1,1)*gmet(1,3)*gmet(2,2)-36*gmet(1,1)*gmet(1,2)*gmet(2,3))&
&       *dgmetds(1,3)+gmet(1,1)*(54*gmet(1,3)**2-18*gmet(1,1)*gmet(3,3))&
&       *dgmetds(2,2)+4*gmet(1,1)*(54*gmet(1,2)*gmet(1,3)-18*gmet(1,1)&
&       *gmet(2,3))*dgmetds(2,3)+gmet(1,1)*(54*gmet(1,2)**2-18*gmet(1,1)&
&       *gmet(2,2))*dgmetds(3,3))/12.d0
       cm(2,12,1,3)=(2*(48*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(6*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+2*gmet(2,2)*(24*gmet(1,3)&
&       *gmet(2,2)+12*gmet(1,2)*gmet(2,3))*dgmetds(1,3)+(-12*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(30*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(2,2)&
&       +4*(6*gmet(1,2)*gmet(1,3)*gmet(2,2)+6*gmet(1,2)**2*gmet(2,3)&
&       +24*gmet(1,1)*gmet(2,2)*gmet(2,3))*dgmetds(2,3)+gmet(2,2)*(12*gmet(1,2)&
&       **2+24*gmet(1,1)*gmet(2,2))*dgmetds(3,3))/4.d0
       cm(3,12,1,3)=(2*gmet(3,3)*(12*gmet(1,3)*gmet(2,3)+24*gmet(1,2)&
&       *gmet(3,3))*dgmetds(1,2)+2*(48*gmet(1,2)*gmet(2,3)*gmet(3,3)&
&       +gmet(1,3)*(6*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,3)&
&       +gmet(3,3)*(12*gmet(1,3)**2+24*gmet(1,1)*gmet(3,3))*dgmetds(2,2)&
&       +4*(6*gmet(1,3)**2*gmet(2,3)+6*gmet(1,2)*gmet(1,3)*gmet(3,3)&
&       +24*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(2,3)+(-12*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(30*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(3,3))&
&       /4.d0
       cm(4,12,1,3)=(2*(6*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(6*gmet(2,3)&
&       **2+24*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+2*(6*gmet(1,3)*gmet(2,2)&
&       *gmet(2,3)+gmet(1,2)*(6*gmet(2,3)**2+24*gmet(2,2)*gmet(3,3)))&
&       *dgmetds(1,3)+(6*gmet(1,3)**2*gmet(2,3)+6*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)+24*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(2,2)+4*(9*gmet(1,3)&
&       **2*gmet(2,2)-6*gmet(1,2)*gmet(1,3)*gmet(2,3)+9*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(9*gmet(2,3)**2+15*gmet(2,2)*gmet(3,3)))*dgmetds(2,3)&
&       +(6*gmet(1,2)*gmet(1,3)*gmet(2,2)+6*gmet(1,2)**2*gmet(2,3)+24*gmet(1,1)&
&       *gmet(2,2)*gmet(2,3))*dgmetds(3,3))/2.d0
       cm(5,12,1,3)=(2*(6*gmet(1,3)**2*gmet(2,3)+48*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)-18*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(1,2)+2*(-12*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)+30*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(-12*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,3)&
&       +(6*gmet(1,3)**3+30*gmet(1,1)*gmet(1,3)*gmet(3,3))*dgmetds(2,2)&
&       +4*(6*gmet(1,1)*gmet(1,3)*gmet(2,3)+gmet(1,2)*(6*gmet(1,3)**2+24*gmet(1,1)&
&       *gmet(3,3)))*dgmetds(2,3)+(6*gmet(1,2)**2*gmet(1,3)-18*gmet(1,1)&
&       *gmet(1,3)*gmet(2,2)+48*gmet(1,1)*gmet(1,2)*gmet(2,3))*dgmetds(3,3))&
&       /4.d0
       cm(6,12,1,3)=(2*(30*gmet(1,3)**2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(-12*gmet(2,3)&
&       **2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+2*(48*gmet(1,2)*gmet(1,3)&
&       *gmet(2,2)+6*gmet(1,2)**2*gmet(2,3)-18*gmet(1,1)*gmet(2,2)*gmet(2,3))&
&       *dgmetds(1,3)+(48*gmet(1,1)*gmet(1,3)*gmet(2,3)+gmet(1,2)*(6*gmet(1,3)&
&       **2-18*gmet(1,1)*gmet(3,3)))*dgmetds(2,2)+4*(6*gmet(1,2)**2*gmet(1,3)&
&       +24*gmet(1,1)*gmet(1,3)*gmet(2,2)+6*gmet(1,1)*gmet(1,2)*gmet(2,3))&
&       *dgmetds(2,3)+(6*gmet(1,2)**3+30*gmet(1,1)*gmet(1,2)*gmet(2,2))&
&       *dgmetds(3,3))/4.d0
       cm(7,12,1,3)=(180*gmet(1,2)*gmet(2,3)**2*dgmetds(2,2)+gmet(2,2)&
&       *(216*gmet(2,3)**2*dgmetds(1,2)-36*gmet(1,2)*gmet(3,3)*dgmetds(2,2)&
&       +gmet(2,3)*(-72*gmet(1,3)*dgmetds(2,2)+432*gmet(1,2)*dgmetds(2,3)))&
&       +gmet(2,2)**2*(-72*gmet(3,3)*dgmetds(1,2)+144*gmet(2,3)*dgmetds(1,3)&
&       -144*gmet(1,3)*dgmetds(2,3)+72*gmet(1,2)*dgmetds(3,3)))/24.d0
       cm(8,12,1,3)=3*gmet(2,3)**3*dgmetds(1,3)+gmet(2,3)**2*(6*gmet(3,3)&
&       *dgmetds(1,2)+6*gmet(1,3)*dgmetds(2,3)+1.5d0*gmet(1,2)*dgmetds(3,3))&
&       +gmet(3,3)*(6*gmet(1,2)*gmet(3,3)*dgmetds(2,2)+gmet(2,2)*(12*gmet(3,3)&
&       *dgmetds(1,2)+24*gmet(1,3)*dgmetds(2,3)-4.5d0*gmet(1,2)*dgmetds(3,3)))&
&       +gmet(2,3)*(gmet(3,3)*(3*gmet(1,3)*dgmetds(2,2)+6*gmet(1,2)*dgmetds(2,3))&
&       +gmet(2,2)*(15*gmet(3,3)*dgmetds(1,3)+12*gmet(1,3)*dgmetds(3,3)))
       cm(9,12,1,3)=3*gmet(2,3)**3*dgmetds(1,2)+gmet(2,3)**2*(6*gmet(2,2)&
&       *dgmetds(1,3)+1.5d0*gmet(1,3)*dgmetds(2,2)+6*gmet(1,2)*dgmetds(2,3))&
&       +gmet(2,3)*(12*gmet(1,2)*gmet(3,3)*dgmetds(2,2)+gmet(2,2)*(15*gmet(3,3)&
&       *dgmetds(1,2)+6*gmet(1,3)*dgmetds(2,3)+3*gmet(1,2)*dgmetds(3,3)))&
&       +gmet(2,2)*(gmet(3,3)*(-4.5d0*gmet(1,3)*dgmetds(2,2)+24*gmet(1,2)&
&       *dgmetds(2,3))+gmet(2,2)*(12*gmet(3,3)*dgmetds(1,3)+6*gmet(1,3)&
&       *dgmetds(3,3)))
       cm(10,12,1,3)=gmet(2,3)*gmet(3,3)*(6*gmet(3,3)*dgmetds(1,2)+18*gmet(1,3)&
&       *dgmetds(2,3)-3*gmet(1,2)*dgmetds(3,3))+gmet(2,3)**2*(9*gmet(3,3)&
&       *dgmetds(1,3)+7.5d0*gmet(1,3)*dgmetds(3,3))+gmet(3,3)*(gmet(3,3)&
&       *(3*gmet(1,3)*dgmetds(2,2)-6*gmet(1,2)*dgmetds(2,3))+gmet(2,2)&
&       *(-3*gmet(3,3)*dgmetds(1,3)-1.5d0*gmet(1,3)*dgmetds(3,3)))
       cm(1,13,1,3)=(180*gmet(1,2)**3*dgmetds(1,3)+gmet(1,1)*gmet(1,2)&
&       *(-216*gmet(2,3)*dgmetds(1,2)-108*gmet(2,2)*dgmetds(1,3)+324*gmet(1,3)&
&       *dgmetds(2,2))+gmet(1,2)**2*(540*gmet(1,3)*dgmetds(1,2)+324*gmet(1,1)&
&       *dgmetds(2,3))-108*gmet(1,1)*(gmet(1,3)*gmet(2,2)*dgmetds(1,2)&
&       +gmet(1,1)*(gmet(2,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(2,3))))&
&       /36.d0
       cm(2,13,1,3)=(72*gmet(1,2)*gmet(2,2)*(gmet(2,3)*dgmetds(1,2)+gmet(2,2)&
&       *dgmetds(1,3))+gmet(1,3)*gmet(2,2)*(144*gmet(2,2)*dgmetds(1,2)&
&       +36*gmet(1,2)*dgmetds(2,2))+144*gmet(1,1)*gmet(2,2)*(gmet(2,3)&
&       *dgmetds(2,2)+gmet(2,2)*dgmetds(2,3))+gmet(1,2)**2*(36*gmet(2,3)&
&       *dgmetds(2,2)+72*gmet(2,2)*dgmetds(2,3)))/12.d0
       cm(3,13,1,3)=(6*(48*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(6*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+2*(-36*gmet(1,3)*gmet(2,2)&
&       *gmet(2,3)+gmet(1,2)*(90*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))&
&       *dgmetds(1,3)+6*(6*gmet(1,3)**2*gmet(2,3)+6*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)+24*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(2,2)+6*(-12*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(30*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(2,3))&
&       /12.d0
       cm(4,13,1,3)=gmet(1,2)*(6*gmet(2,3)**2*dgmetds(1,2)+24*gmet(2,2)&
&       *gmet(3,3)*dgmetds(1,2)+18*gmet(2,2)*gmet(2,3)*dgmetds(1,3))&
&       +9*gmet(1,3)**2*gmet(2,2)*dgmetds(2,2)+gmet(1,2)**2*(9*gmet(3,3)&
&       *dgmetds(2,2)+6*gmet(2,3)*dgmetds(2,3))+gmet(1,1)*(9*gmet(2,3)&
&       **2*dgmetds(2,2)+15*gmet(2,2)*gmet(3,3)*dgmetds(2,2)+24*gmet(2,2)&
&       *gmet(2,3)*dgmetds(2,3))+gmet(1,3)*(-6*gmet(2,2)**2*dgmetds(1,3)&
&       -6*gmet(1,2)*gmet(2,3)*dgmetds(2,2)+6*gmet(2,2)*(gmet(2,3)*dgmetds(1,2)&
&       +gmet(1,2)*dgmetds(2,3)))
       cm(5,13,1,3)=(6*(-12*gmet(1,3)**2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)+30*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(-12*gmet(2,3)&
&       **2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+2*(-36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,2)+90*gmet(1,2)**2*gmet(2,3)-18*gmet(1,1)*gmet(2,2)*gmet(2,3))&
&       *dgmetds(1,3)+6*(6*gmet(1,1)*gmet(1,3)*gmet(2,3)+gmet(1,2)*(6*gmet(1,3)&
&       **2+24*gmet(1,1)*gmet(3,3)))*dgmetds(2,2)+6*(6*gmet(1,2)**2*gmet(1,3)&
&       -18*gmet(1,1)*gmet(1,3)*gmet(2,2)+48*gmet(1,1)*gmet(1,2)*gmet(2,3))&
&       *dgmetds(2,3))/12.d0
       cm(6,13,1,3)=gmet(1,2)**2*(3*gmet(2,3)*dgmetds(1,2)+9*gmet(2,2)&
&       *dgmetds(1,3)+3*gmet(1,3)*dgmetds(2,2))+gmet(1,1)*gmet(2,2)*(-9*gmet(2,3)&
&       *dgmetds(1,2)-3*gmet(2,2)*dgmetds(1,3)+12*gmet(1,3)*dgmetds(2,2))&
&       +3*gmet(1,2)**3*dgmetds(2,3)+gmet(1,2)*(24*gmet(1,3)*gmet(2,2)&
&       *dgmetds(1,2)+gmet(1,1)*(3*gmet(2,3)*dgmetds(2,2)+15*gmet(2,2)&
&       *dgmetds(2,3)))
       cm(7,13,1,3)=gmet(2,2)*(2*gmet(2,2)**2*dgmetds(1,3)+9*gmet(1,2)&
&       *gmet(2,3)*dgmetds(2,2)+gmet(2,2)*(6*gmet(2,3)*dgmetds(1,2)-3*gmet(1,3)&
&       *dgmetds(2,2)+6*gmet(1,2)*dgmetds(2,3)))
       cm(8,13,1,3)=3*gmet(2,3)**3*dgmetds(1,2)+gmet(2,2)*gmet(3,3)*(-3*gmet(2,2)&
&       *dgmetds(1,3)+12*gmet(1,3)*dgmetds(2,2)-9*gmet(1,2)*dgmetds(2,3))&
&       +gmet(2,3)**2*(9*gmet(2,2)*dgmetds(1,3)+3*(gmet(1,3)*dgmetds(2,2)&
&       +gmet(1,2)*dgmetds(2,3)))+gmet(2,3)*(3*gmet(1,2)*gmet(3,3)*dgmetds(2,2)&
&       +gmet(2,2)*(15*gmet(3,3)*dgmetds(1,2)+24*gmet(1,3)*dgmetds(2,3)))
       cm(9,13,1,3)=3*gmet(1,2)*gmet(2,3)**2*dgmetds(2,2)+gmet(2,2)**2*(12*gmet(3,3)&
&       *dgmetds(1,2)+6*gmet(2,3)*dgmetds(1,3)+12*gmet(1,3)*dgmetds(2,3))&
&       +gmet(2,2)*(6*gmet(2,3)**2*dgmetds(1,2)+12*gmet(1,2)*gmet(3,3)&
&       *dgmetds(2,2)+gmet(2,3)*(3*gmet(1,3)*dgmetds(2,2)+6*gmet(1,2)&
&       *dgmetds(2,3)))
       cm(10,13,1,3)=5*gmet(2,3)**3*dgmetds(1,3)+gmet(2,3)*gmet(3,3)&
&       *(-3*gmet(2,2)*dgmetds(1,3)+9*gmet(1,3)*dgmetds(2,2)-6*gmet(1,2)&
&       *dgmetds(2,3))+gmet(2,3)**2*(9*gmet(3,3)*dgmetds(1,2)+15*gmet(1,3)&
&       *dgmetds(2,3))-3*gmet(3,3)*(gmet(1,2)*gmet(3,3)*dgmetds(2,2)&
&       +gmet(2,2)*(gmet(3,3)*dgmetds(1,2)+gmet(1,3)*dgmetds(2,3)))
       cm(1,14,1,3)=(180*gmet(1,3)**3*dgmetds(1,2)+gmet(1,3)**2*(540*gmet(1,2)&
&       *dgmetds(1,3)+324*gmet(1,1)*dgmetds(2,3))+gmet(1,1)*gmet(1,3)&
&       *(-108*gmet(3,3)*dgmetds(1,2)-216*gmet(2,3)*dgmetds(1,3)+324*gmet(1,2)&
&       *dgmetds(3,3))-108*gmet(1,1)*(gmet(1,2)*gmet(3,3)*dgmetds(1,3)&
&       +gmet(1,1)*(gmet(3,3)*dgmetds(2,3)+gmet(2,3)*dgmetds(3,3))))&
&       /36.d0
       cm(2,14,1,3)=(2*(-36*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+6*(48*gmet(1,3)*gmet(2,2)&
&       *gmet(2,3)+gmet(1,2)*(6*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))&
&       *dgmetds(1,3)+6*(-12*gmet(1,3)**2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(30*gmet(2,3)&
&       **2-6*gmet(2,2)*gmet(3,3)))*dgmetds(2,3)+6*(6*gmet(1,2)*gmet(1,3)&
&       *gmet(2,2)+6*gmet(1,2)**2*gmet(2,3)+24*gmet(1,1)*gmet(2,2)*gmet(2,3))&
&       *dgmetds(3,3))/12.d0
       cm(3,14,1,3)=(gmet(1,3)*gmet(3,3)*(72*gmet(3,3)*dgmetds(1,2)+72*gmet(2,3)&
&       *dgmetds(1,3)+36*gmet(1,2)*dgmetds(3,3))+gmet(1,3)**2*(72*gmet(3,3)&
&       *dgmetds(2,3)+36*gmet(2,3)*dgmetds(3,3))+144*gmet(3,3)*(gmet(1,2)&
&       *gmet(3,3)*dgmetds(1,3)+gmet(1,1)*(gmet(3,3)*dgmetds(2,3)+gmet(2,3)&
&       *dgmetds(3,3))))/12.d0
       cm(4,14,1,3)=(2*gmet(3,3)*(54*gmet(1,3)*gmet(2,3)-18*gmet(1,2)&
&       *gmet(3,3))*dgmetds(1,2)+6*(6*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)&
&       *(6*gmet(2,3)**2+24*gmet(2,2)*gmet(3,3)))*dgmetds(1,3)+6*(6*gmet(1,3)&
&       **2*gmet(2,3)+6*gmet(1,2)*gmet(1,3)*gmet(3,3)+24*gmet(1,1)*gmet(2,3)&
&       *gmet(3,3))*dgmetds(2,3)+6*(9*gmet(1,3)**2*gmet(2,2)-6*gmet(1,2)&
&       *gmet(1,3)*gmet(2,3)+9*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(9*gmet(2,3)&
&       **2+15*gmet(2,2)*gmet(3,3)))*dgmetds(3,3))/6.d0
       cm(5,14,1,3)=(36*gmet(1,3)**3*dgmetds(2,3)+gmet(1,1)*gmet(3,3)&
&       *(-36*gmet(3,3)*dgmetds(1,2)-108*gmet(2,3)*dgmetds(1,3)+144*gmet(1,2)&
&       *dgmetds(3,3))+gmet(1,3)**2*(108*gmet(3,3)*dgmetds(1,2)+36*(gmet(2,3)&
&       *dgmetds(1,3)+gmet(1,2)*dgmetds(3,3)))+gmet(1,3)*(288*gmet(1,2)&
&       *gmet(3,3)*dgmetds(1,3)+gmet(1,1)*(180*gmet(3,3)*dgmetds(2,3)&
&       +36*gmet(2,3)*dgmetds(3,3))))/12.d0
       cm(6,14,1,3)=(2*(90*gmet(1,3)**2*gmet(2,3)-36*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)-18*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(1,2)+6*(30*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(-12*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,3)&
&       +6*(48*gmet(1,1)*gmet(1,3)*gmet(2,3)+gmet(1,2)*(6*gmet(1,3)**2-18*gmet(1,1)&
&       *gmet(3,3)))*dgmetds(2,3)+6*(6*gmet(1,2)**2*gmet(1,3)+24*gmet(1,1)&
&       *gmet(1,3)*gmet(2,2)+6*gmet(1,1)*gmet(1,2)*gmet(2,3))*dgmetds(3,3))&
&       /12.d0
       cm(7,14,1,3)=(180*gmet(2,3)**3*dgmetds(1,2)+gmet(2,3)**2*(324*gmet(2,2)&
&       *dgmetds(1,3)+540*gmet(1,2)*dgmetds(2,3))+gmet(2,2)*gmet(2,3)&
&       *(-108*gmet(3,3)*dgmetds(1,2)-216*gmet(1,3)*dgmetds(2,3)+324*gmet(1,2)&
&       *dgmetds(3,3))-108*gmet(2,2)*(gmet(1,2)*gmet(3,3)*dgmetds(2,3)&
&       +gmet(2,2)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))))&
&       /36.d0
       cm(8,14,1,3)=(gmet(2,3)*gmet(3,3)*(72*gmet(3,3)*dgmetds(1,2)+72*gmet(1,3)&
&       *dgmetds(2,3)+36*gmet(1,2)*dgmetds(3,3))+gmet(2,3)**2*(72*gmet(3,3)&
&       *dgmetds(1,3)+36*gmet(1,3)*dgmetds(3,3))+144*gmet(3,3)*(gmet(1,2)&
&       *gmet(3,3)*dgmetds(2,3)+gmet(2,2)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)&
&       *dgmetds(3,3))))/12.d0
       cm(9,14,1,3)=(36*gmet(2,3)**3*dgmetds(1,3)+gmet(2,2)*gmet(3,3)&
&       *(-36*gmet(3,3)*dgmetds(1,2)-108*gmet(1,3)*dgmetds(2,3)+144*gmet(1,2)&
&       *dgmetds(3,3))+gmet(2,3)**2*(108*gmet(3,3)*dgmetds(1,2)+36*(gmet(1,3)&
&       *dgmetds(2,3)+gmet(1,2)*dgmetds(3,3)))+gmet(2,3)*(288*gmet(1,2)&
&       *gmet(3,3)*dgmetds(2,3)+gmet(2,2)*(180*gmet(3,3)*dgmetds(1,3)&
&       +36*gmet(1,3)*dgmetds(3,3))))/12.d0
       cm(10,14,1,3)=gmet(3,3)*(2*gmet(3,3)**2*dgmetds(1,2)+9*gmet(1,3)&
&       *gmet(2,3)*dgmetds(3,3)+gmet(3,3)*(6*gmet(2,3)*dgmetds(1,3)+6*gmet(1,3)&
&       *dgmetds(2,3)-3*gmet(1,2)*dgmetds(3,3)))
       cm(1,15,1,3)=5*gmet(1,3)**3*dgmetds(1,3)-3*gmet(1,1)*gmet(1,3)&
&       *gmet(3,3)*dgmetds(1,3)+4.5d0*gmet(1,1)*gmet(1,3)**2*dgmetds(3,3)&
&       -1.5d0*gmet(1,1)**2*gmet(3,3)*dgmetds(3,3)
       cm(2,15,1,3)=(8*(-36*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,3)+12*(-12*gmet(1,3)**2*gmet(2,2)&
&       +36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)&
&       *(30*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(3,3))/48.d0
       cm(3,15,1,3)=gmet(3,3)*(6*gmet(1,3)*gmet(3,3)*dgmetds(1,3)+3*gmet(1,3)&
&       **2*dgmetds(3,3)+6*gmet(1,1)*gmet(3,3)*dgmetds(3,3))
       cm(4,15,1,3)=3*gmet(1,3)**2*gmet(2,3)*dgmetds(3,3)+gmet(1,3)*gmet(3,3)&
&       *(18*gmet(2,3)*dgmetds(1,3)+3*gmet(1,2)*dgmetds(3,3))+gmet(3,3)&
&       *(-6*gmet(1,2)*gmet(3,3)*dgmetds(1,3)+12*gmet(1,1)*gmet(2,3)&
&       *dgmetds(3,3))
       cm(5,15,1,3)=9*gmet(1,3)**2*gmet(3,3)*dgmetds(1,3)-3*gmet(1,1)&
&       *gmet(3,3)**2*dgmetds(1,3)+1.5d0*gmet(1,3)**3*dgmetds(3,3)+7.5d0*gmet(1,1)&
&       *gmet(1,3)*gmet(3,3)*dgmetds(3,3)
       cm(6,15,1,3)=(8*(90*gmet(1,3)**2*gmet(2,3)-36*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)-18*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(1,3)+12*(48*gmet(1,1)&
&       *gmet(1,3)*gmet(2,3)+gmet(1,2)*(6*gmet(1,3)**2-18*gmet(1,1)*gmet(3,3)))&
&       *dgmetds(3,3))/48.d0
       cm(7,15,1,3)=5*gmet(2,3)**3*dgmetds(1,3)+7.5d0*gmet(1,2)*gmet(2,3)&
&       **2*dgmetds(3,3)-1.5d0*gmet(1,2)*gmet(2,2)*gmet(3,3)*dgmetds(3,3)&
&       -3*gmet(2,2)*gmet(2,3)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(8,15,1,3)=gmet(3,3)*(6*gmet(1,2)*gmet(3,3)*dgmetds(3,3)+gmet(2,3)&
&       *(6*gmet(3,3)*dgmetds(1,3)+3*gmet(1,3)*dgmetds(3,3)))
       cm(9,15,1,3)=12*gmet(1,2)*gmet(2,3)*gmet(3,3)*dgmetds(3,3)+gmet(2,2)&
&       *gmet(3,3)*(-3*gmet(3,3)*dgmetds(1,3)-4.5d0*gmet(1,3)*dgmetds(3,3))&
&       +gmet(2,3)**2*(9*gmet(3,3)*dgmetds(1,3)+1.5d0*gmet(1,3)*dgmetds(3,3))
       cm(10,15,1,3)=gmet(3,3)**2*(2*gmet(3,3)*dgmetds(1,3)+3*gmet(1,3)&
&       *dgmetds(3,3))
       cm(1,16,1,3)=gmet(1,2)*(2.5d0*gmet(1,2)**2-1.5d0*gmet(1,1)*gmet(2,2))&
&       *dgmetds(2,2)
       cm(2,16,1,3)=3*gmet(1,2)*gmet(2,2)**2*dgmetds(2,2)
       cm(3,16,1,3)=((-36*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(2,2))/12.d0
       cm(4,16,1,3)=(gmet(2,2)*(-18*gmet(1,3)*gmet(2,2)+54*gmet(1,2)&
&       *gmet(2,3))*dgmetds(2,2))/6.d0
       cm(5,16,1,3)=((-36*gmet(1,2)*gmet(1,3)*gmet(2,2)+90*gmet(1,2)&
&       **2*gmet(2,3)-18*gmet(1,1)*gmet(2,2)*gmet(2,3))*dgmetds(2,2))&
&       /12.d0
       cm(6,16,1,3)=gmet(2,2)*(4.5d0*gmet(1,2)**2-1.5d0*gmet(1,1)*gmet(2,2))&
&       *dgmetds(2,2)
       cm(7,16,1,3)=gmet(2,2)**3*dgmetds(2,2)
       cm(8,16,1,3)=gmet(2,2)*(4.5d0*gmet(2,3)**2-1.5d0*gmet(2,2)*gmet(3,3))&
&       *dgmetds(2,2)
       cm(9,16,1,3)=3*gmet(2,2)**2*gmet(2,3)*dgmetds(2,2)
       cm(10,16,1,3)=gmet(2,3)*(2.5d0*gmet(2,3)**2-1.5d0*gmet(2,2)*gmet(3,3))&
&       *dgmetds(2,2)
       cm(1,17,1,3)=(1080*gmet(1,2)**2*gmet(1,3)*dgmetds(2,3)-216*gmet(1,1)&
&       *gmet(1,3)*(gmet(2,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(2,3))+180*gmet(1,2)&
&       **3*dgmetds(3,3)+gmet(1,2)*(540*gmet(1,3)**2*dgmetds(2,2)+gmet(1,1)&
&       *(-108*gmet(3,3)*dgmetds(2,2)-432*gmet(2,3)*dgmetds(2,3)-108*gmet(2,2)&
&       *dgmetds(3,3))))/72.d0
       cm(2,17,1,3)=(288*gmet(1,3)*gmet(2,2)*(gmet(2,3)*dgmetds(2,2)&
&       +gmet(2,2)*dgmetds(2,3))+gmet(1,2)*(36*gmet(2,3)**2*dgmetds(2,2)&
&       +144*gmet(2,2)*gmet(2,3)*dgmetds(2,3)+gmet(2,2)*(-108*gmet(3,3)&
&       *dgmetds(2,2)+72*gmet(2,2)*dgmetds(3,3))))/24.d0
       cm(3,17,1,3)=gmet(1,3)*(3*gmet(2,3)**2*dgmetds(2,3)-9*gmet(2,2)&
&       *gmet(3,3)*dgmetds(2,3)+gmet(2,3)*(3*gmet(3,3)*dgmetds(2,2)-3*gmet(2,2)&
&       *dgmetds(3,3)))+gmet(1,2)*(6*gmet(3,3)**2*dgmetds(2,2)+7.5d0*gmet(2,3)&
&       **2*dgmetds(3,3)+gmet(3,3)*(24*gmet(2,3)*dgmetds(2,3)-1.5d0*gmet(2,2)&
&       *dgmetds(3,3)))
       cm(4,17,1,3)=gmet(1,3)*(3*gmet(2,3)**2*dgmetds(2,2)+6*gmet(2,2)&
&       *gmet(2,3)*dgmetds(2,3)+gmet(2,2)*(12*gmet(3,3)*dgmetds(2,2)&
&       -3*gmet(2,2)*dgmetds(3,3)))+gmet(1,2)*(6*gmet(2,3)**2*dgmetds(2,3)&
&       +24*gmet(2,2)*gmet(3,3)*dgmetds(2,3)+gmet(2,3)*(3*gmet(3,3)*dgmetds(2,2)&
&       +9*gmet(2,2)*dgmetds(3,3)))
       cm(5,17,1,3)=(6*(6*gmet(1,3)**2*gmet(2,3)+48*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)-18*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(2,2)+12*(-12*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)+30*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(-12*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(2,3)&
&       +2*(-36*gmet(1,2)*gmet(1,3)*gmet(2,2)+90*gmet(1,2)**2*gmet(2,3)&
&       -18*gmet(1,1)*gmet(2,2)*gmet(2,3))*dgmetds(3,3))/24.d0
       cm(6,17,1,3)=7.5d0*gmet(1,3)**2*gmet(2,2)*dgmetds(2,2)+gmet(1,2)&
&       *gmet(1,3)*(9*gmet(2,3)*dgmetds(2,2)+24*gmet(2,2)*dgmetds(2,3))&
&       +gmet(1,2)**2*(-3*gmet(3,3)*dgmetds(2,2)+3*gmet(2,3)*dgmetds(2,3)&
&       +4.5d0*gmet(2,2)*dgmetds(3,3))+gmet(1,1)*(-3*gmet(2,3)**2*dgmetds(2,2)&
&       -9*gmet(2,2)*gmet(2,3)*dgmetds(2,3)-1.5d0*gmet(2,2)*(gmet(3,3)&
&       *dgmetds(2,2)+gmet(2,2)*dgmetds(3,3)))
       cm(7,17,1,3)=gmet(2,2)*(4.5d0*gmet(2,3)**2*dgmetds(2,2)+6*gmet(2,2)&
&       *gmet(2,3)*dgmetds(2,3)+gmet(2,2)*(-1.5d0*gmet(3,3)*dgmetds(2,2)&
&       +gmet(2,2)*dgmetds(3,3)))
       cm(8,17,1,3)=3*gmet(2,3)**3*dgmetds(2,3)+15*gmet(2,2)*gmet(2,3)&
&       *gmet(3,3)*dgmetds(2,3)+gmet(2,2)*gmet(3,3)*(6*gmet(3,3)*dgmetds(2,2)&
&       -1.5d0*gmet(2,2)*dgmetds(3,3))+gmet(2,3)**2*(3*gmet(3,3)*dgmetds(2,2)&
&       +4.5d0*gmet(2,2)*dgmetds(3,3))
       cm(9,17,1,3)=1.5d0*gmet(2,3)**3*dgmetds(2,2)+6*gmet(2,2)*gmet(2,3)&
&       **2*dgmetds(2,3)+12*gmet(2,2)**2*gmet(3,3)*dgmetds(2,3)+gmet(2,2)&
&       *gmet(2,3)*(7.5d0*gmet(3,3)*dgmetds(2,2)+3*gmet(2,2)*dgmetds(3,3))
       cm(10,17,1,3)=9*gmet(2,3)**2*gmet(3,3)*dgmetds(2,3)-3*gmet(2,2)&
&       *gmet(3,3)**2*dgmetds(2,3)+2.5d0*gmet(2,3)**3*dgmetds(3,3)+gmet(2,3)&
&       *gmet(3,3)*(3*gmet(3,3)*dgmetds(2,2)-1.5d0*gmet(2,2)*dgmetds(3,3))
       cm(1,18,1,3)=7.5d0*gmet(1,2)**2*gmet(1,3)*dgmetds(2,2)-1.5d0*gmet(1,1)&
&       *gmet(1,3)*gmet(2,2)*dgmetds(2,2)+5*gmet(1,2)**3*dgmetds(2,3)&
&       -3*gmet(1,1)*gmet(1,2)*(gmet(2,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(2,3))
       cm(2,18,1,3)=gmet(2,2)*(6*gmet(1,3)*gmet(2,2)*dgmetds(2,2)+gmet(1,2)&
&       *(3*gmet(2,3)*dgmetds(2,2)+6*gmet(2,2)*dgmetds(2,3)))
       cm(3,18,1,3)=gmet(1,3)*(1.5d0*gmet(2,3)**2*dgmetds(2,2)-4.5d0*gmet(2,2)&
&       *gmet(3,3)*dgmetds(2,2)-6*gmet(2,2)*gmet(2,3)*dgmetds(2,3))+gmet(1,2)&
&       *(12*gmet(2,3)*gmet(3,3)*dgmetds(2,2)+15*gmet(2,3)**2*dgmetds(2,3)&
&       -3*gmet(2,2)*gmet(3,3)*dgmetds(2,3))
       cm(4,18,1,3)=gmet(1,3)*gmet(2,2)*(3*gmet(2,3)*dgmetds(2,2)-6*gmet(2,2)&
&       *dgmetds(2,3))+gmet(1,2)*(3*gmet(2,3)**2*dgmetds(2,2)+12*gmet(2,2)&
&       *gmet(3,3)*dgmetds(2,2)+18*gmet(2,2)*gmet(2,3)*dgmetds(2,3))
       cm(5,18,1,3)=(12*(-12*gmet(1,3)**2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)+30*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(-12*gmet(2,3)&
&       **2-6*gmet(2,2)*gmet(3,3)))*dgmetds(2,2)+8*(-36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,2)+90*gmet(1,2)**2*gmet(2,3)-18*gmet(1,1)*gmet(2,2)*gmet(2,3))&
&       *dgmetds(2,3))/48.d0
       cm(6,18,1,3)=12*gmet(1,2)*gmet(1,3)*gmet(2,2)*dgmetds(2,2)+gmet(1,1)&
&       *gmet(2,2)*(-4.5d0*gmet(2,3)*dgmetds(2,2)-3*gmet(2,2)*dgmetds(2,3))&
&       +gmet(1,2)**2*(1.5d0*gmet(2,3)*dgmetds(2,2)+9*gmet(2,2)*dgmetds(2,3))
       cm(7,18,1,3)=gmet(2,2)**2*(3*gmet(2,3)*dgmetds(2,2)+2*gmet(2,2)&
&       *dgmetds(2,3))
       cm(8,18,1,3)=1.5d0*gmet(2,3)**3*dgmetds(2,2)+7.5d0*gmet(2,2)*gmet(2,3)&
&       *gmet(3,3)*dgmetds(2,2)+9*gmet(2,2)*gmet(2,3)**2*dgmetds(2,3)&
&       -3*gmet(2,2)**2*gmet(3,3)*dgmetds(2,3)
       cm(9,18,1,3)=gmet(2,2)*(3*gmet(2,3)**2*dgmetds(2,2)+6*gmet(2,2)&
&       *gmet(3,3)*dgmetds(2,2)+6*gmet(2,2)*gmet(2,3)*dgmetds(2,3))
       cm(10,18,1,3)=4.5d0*gmet(2,3)**2*gmet(3,3)*dgmetds(2,2)-1.5d0*gmet(2,2)&
&       *gmet(3,3)**2*dgmetds(2,2)+5*gmet(2,3)**3*dgmetds(2,3)-3*gmet(2,2)&
&       *gmet(2,3)*gmet(3,3)*dgmetds(2,3)
       cm(1,19,1,3)=(180*gmet(1,3)**3*dgmetds(2,2)+1080*gmet(1,2)*gmet(1,3)&
&       **2*dgmetds(2,3)-216*gmet(1,1)*gmet(1,2)*(gmet(3,3)*dgmetds(2,3)&
&       +gmet(2,3)*dgmetds(3,3))+gmet(1,3)*(540*gmet(1,2)**2*dgmetds(3,3)&
&       +gmet(1,1)*(-108*gmet(3,3)*dgmetds(2,2)-432*gmet(2,3)*dgmetds(2,3)&
&       -108*gmet(2,2)*dgmetds(3,3))))/72.d0
       cm(2,19,1,3)=(2*(-36*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(2,2)+12*(48*gmet(1,3)*gmet(2,2)&
&       *gmet(2,3)+gmet(1,2)*(6*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))&
&       *dgmetds(2,3)+6*gmet(2,2)*(24*gmet(1,3)*gmet(2,2)+12*gmet(1,2)&
&       *gmet(2,3))*dgmetds(3,3))/24.0d0
       cm(3,19,1,3)=12*gmet(1,2)*gmet(3,3)*(gmet(3,3)*dgmetds(2,3)+gmet(2,3)&
&       *dgmetds(3,3))+gmet(1,3)*(3*gmet(3,3)**2*dgmetds(2,2)+1.5d0*gmet(2,3)&
&       **2*dgmetds(3,3)+gmet(3,3)*(6*gmet(2,3)*dgmetds(2,3)-4.5d0*gmet(2,2)&
&       *dgmetds(3,3)))
       cm(4,19,1,3)=gmet(1,3)*(6*gmet(2,3)**2*dgmetds(2,3)+24*gmet(2,2)&
&       *gmet(3,3)*dgmetds(2,3)+gmet(2,3)*(9*gmet(3,3)*dgmetds(2,2)+3*gmet(2,2)&
&       *dgmetds(3,3)))+gmet(1,2)*(-3*gmet(3,3)**2*dgmetds(2,2)+3*gmet(2,3)&
&       **2*dgmetds(3,3)+gmet(3,3)*(6*gmet(2,3)*dgmetds(2,3)+12*gmet(2,2)&
&       *dgmetds(3,3)))
       cm(5,19,1,3)=(2*gmet(3,3)*(54*gmet(1,3)**2-18*gmet(1,1)*gmet(3,3))&
&       *dgmetds(2,2)+12*(6*gmet(1,3)**2*gmet(2,3)+48*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)-18*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(2,3)+6*(-12*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)+30*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(-12*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(3,3))&
&       /24.d0
       cm(6,19,1,3)=gmet(1,3)**2*(7.5d0*gmet(2,3)*dgmetds(2,2)+15*gmet(2,2)&
&       *dgmetds(2,3))+gmet(1,2)*gmet(1,3)*(-3*gmet(3,3)*dgmetds(2,2)&
&       +18*gmet(2,3)*dgmetds(2,3)+12*gmet(2,2)*dgmetds(3,3))+gmet(1,2)&
&       **2*(-6*gmet(3,3)*dgmetds(2,3)+1.5d0*gmet(2,3)*dgmetds(3,3))&
&       +gmet(1,1)*(-1.5d0*gmet(2,3)*gmet(3,3)*dgmetds(2,2)-6*gmet(2,3)&
&       **2*dgmetds(2,3)-3*gmet(2,2)*gmet(3,3)*dgmetds(2,3)-4.5d0*gmet(2,2)&
&       *gmet(2,3)*dgmetds(3,3))
       cm(7,19,1,3)=2.5d0*gmet(2,3)**3*dgmetds(2,2)+9*gmet(2,2)*gmet(2,3)&
&       **2*dgmetds(2,3)-3*gmet(2,2)**2*gmet(3,3)*dgmetds(2,3)+gmet(2,2)&
&       *gmet(2,3)*(-1.5d0*gmet(3,3)*dgmetds(2,2)+3*gmet(2,2)*dgmetds(3,3))
       cm(8,19,1,3)=6*gmet(2,3)**2*gmet(3,3)*dgmetds(2,3)+12*gmet(2,2)&
&       *gmet(3,3)**2*dgmetds(2,3)+1.5d0*gmet(2,3)**3*dgmetds(3,3)+gmet(2,3)&
&       *gmet(3,3)*(3*gmet(3,3)*dgmetds(2,2)+7.5d0*gmet(2,2)*dgmetds(3,3))
       cm(9,19,1,3)=3*gmet(2,3)**3*dgmetds(2,3)+15*gmet(2,2)*gmet(2,3)&
&       *gmet(3,3)*dgmetds(2,3)+gmet(2,3)**2*(4.5d0*gmet(3,3)*dgmetds(2,2)&
&       +3*gmet(2,2)*dgmetds(3,3))+gmet(2,2)*gmet(3,3)*(-1.5d0*gmet(3,3)&
&       *dgmetds(2,2)+6*gmet(2,2)*dgmetds(3,3))
       cm(10,19,1,3)=gmet(3,3)*(1*gmet(3,3)**2*dgmetds(2,2)+4.5d0*gmet(2,3)&
&       **2*dgmetds(3,3)+gmet(3,3)*(6*gmet(2,3)*dgmetds(2,3)-1.5d0*gmet(2,2)&
&       *dgmetds(3,3)))
       cm(1,20,1,3)=5*gmet(1,3)**3*dgmetds(2,3)+7.5d0*gmet(1,2)*gmet(1,3)&
&       **2*dgmetds(3,3)-1.5d0*gmet(1,1)*gmet(1,2)*gmet(3,3)*dgmetds(3,3)&
&       -3*gmet(1,1)*gmet(1,3)*(gmet(3,3)*dgmetds(2,3)+gmet(2,3)*dgmetds(3,3))
       cm(2,20,1,3)=(8*(-36*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(2,3)+12*(48*gmet(1,3)*gmet(2,2)&
&       *gmet(2,3)+gmet(1,2)*(6*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))&
&       *dgmetds(3,3))/48.d0
       cm(3,20,1,3)=gmet(3,3)*(6*gmet(1,2)*gmet(3,3)*dgmetds(3,3)+gmet(1,3)&
&       *(6*gmet(3,3)*dgmetds(2,3)+3*gmet(2,3)*dgmetds(3,3)))
       cm(4,20,1,3)=gmet(1,2)*gmet(3,3)*(-6*gmet(3,3)*dgmetds(2,3)+3*gmet(2,3)&
&       *dgmetds(3,3))+gmet(1,3)*(18*gmet(2,3)*gmet(3,3)*dgmetds(2,3)&
&       +3*gmet(2,3)**2*dgmetds(3,3)+12*gmet(2,2)*gmet(3,3)*dgmetds(3,3))
       cm(5,20,1,3)=12*gmet(1,2)*gmet(1,3)*gmet(3,3)*dgmetds(3,3)+gmet(1,1)&
&       *gmet(3,3)*(-3*gmet(3,3)*dgmetds(2,3)-4.5d0*gmet(2,3)*dgmetds(3,3))&
&       +gmet(1,3)**2*(9*gmet(3,3)*dgmetds(2,3)+1.5d0*gmet(2,3)*dgmetds(3,3))
       cm(6,20,1,3)=(8*(90*gmet(1,3)**2*gmet(2,3)-36*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)-18*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(2,3)+12*(30*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(-12*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(3,3))&
&       /48.d0
       cm(7,20,1,3)=5*gmet(2,3)**3*dgmetds(2,3)-3*gmet(2,2)*gmet(2,3)&
&       *gmet(3,3)*dgmetds(2,3)+4.5d0*gmet(2,2)*gmet(2,3)**2*dgmetds(3,3)&
&       -1.5d0*gmet(2,2)**2*gmet(3,3)*dgmetds(3,3)
       cm(8,20,1,3)=gmet(3,3)*(6*gmet(2,3)*gmet(3,3)*dgmetds(2,3)+3*gmet(2,3)&
&       **2*dgmetds(3,3)+6*gmet(2,2)*gmet(3,3)*dgmetds(3,3))
       cm(9,20,1,3)=9*gmet(2,3)**2*gmet(3,3)*dgmetds(2,3)-3*gmet(2,2)&
&       *gmet(3,3)**2*dgmetds(2,3)+1.5d0*gmet(2,3)**3*dgmetds(3,3)+7.5d0*gmet(2,2)&
&       *gmet(2,3)*gmet(3,3)*dgmetds(3,3)
       cm(10,20,1,3)=gmet(3,3)**2*(2*gmet(3,3)*dgmetds(2,3)+3*gmet(2,3)&
&       *dgmetds(3,3))
       cm(1,21,1,3)=gmet(1,3)*(2.5d0*gmet(1,3)**2-1.5d0*gmet(1,1)*gmet(3,3))&
&       *dgmetds(3,3)
       cm(2,21,1,3)=((-36*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(3,3))/12.d0
       cm(3,21,1,3)=3*gmet(1,3)*gmet(3,3)**2*dgmetds(3,3)
       cm(4,21,1,3)=(gmet(3,3)*(54*gmet(1,3)*gmet(2,3)-18*gmet(1,2)*gmet(3,3))&
&       *dgmetds(3,3))/6.d0
       cm(5,21,1,3)=gmet(3,3)*(4.5d0*gmet(1,3)**2-1.5d0*gmet(1,1)*gmet(3,3))&
&       *dgmetds(3,3)
       cm(6,21,1,3)=((90*gmet(1,3)**2*gmet(2,3)-36*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)-18*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(3,3))/12.d0
       cm(7,21,1,3)=gmet(2,3)*(2.5d0*gmet(2,3)**2-1.5d0*gmet(2,2)*gmet(3,3))&
&       *dgmetds(3,3)
       cm(8,21,1,3)=3*gmet(2,3)*gmet(3,3)**2*dgmetds(3,3)
       cm(9,21,1,3)=gmet(3,3)*(4.5d0*gmet(2,3)**2-1.5d0*gmet(2,2)*gmet(3,3))&
&       *dgmetds(3,3)
       cm(10,21,1,3)=gmet(3,3)**3*dgmetds(3,3)
     elseif(iterm==2)then
       cm(1,1,2,3)=6*gmet(1,1)**2*dgmetds(1,1)
       cm(2,1,2,3)=9*gmet(1,2)**2*dgmetds(1,1)+18*gmet(1,1)*gmet(1,2)&
&       *dgmetds(1,2)+gmet(1,1)*(-6*gmet(2,2)*dgmetds(1,1)-3*gmet(1,1)&
&       *dgmetds(2,2))
       cm(3,1,2,3)=9*gmet(1,3)**2*dgmetds(1,1)+18*gmet(1,1)*gmet(1,3)&
&       *dgmetds(1,3)+gmet(1,1)*(-6*gmet(3,3)*dgmetds(1,1)-3*gmet(1,1)&
&       *dgmetds(3,3))
       cm(4,1,2,3)=18*gmet(1,2)*(gmet(1,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,3))&
&       +gmet(1,1)*(-12*gmet(2,3)*dgmetds(1,1)+18*gmet(1,3)*dgmetds(1,2)&
&       -6*gmet(1,1)*dgmetds(2,3))
       cm(5,1,2,3)=gmet(1,1)*(12*gmet(1,3)*dgmetds(1,1)+6*gmet(1,1)*dgmetds(1,3))
       cm(6,1,2,3)=gmet(1,1)*(12*gmet(1,2)*dgmetds(1,1)+6*gmet(1,1)*dgmetds(1,2))
       cm(7,1,2,3)=15*gmet(1,2)**2*dgmetds(1,2)-3*gmet(1,1)*gmet(2,2)&
&       *dgmetds(1,2)-3*gmet(1,2)*(gmet(2,2)*dgmetds(1,1)+gmet(1,1)*dgmetds(2,2))
       cm(8,1,2,3)=15*gmet(1,3)**2*dgmetds(1,2)+gmet(1,1)*(-3*gmet(3,3)&
&       *dgmetds(1,2)-6*gmet(2,3)*dgmetds(1,3))+gmet(1,3)*(-6*gmet(2,3)&
&       *dgmetds(1,1)+30*gmet(1,2)*dgmetds(1,3)-6*gmet(1,1)*dgmetds(2,3))&
&       -3*gmet(1,2)*(gmet(3,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(3,3))
       cm(9,1,2,3)=15*gmet(1,2)**2*dgmetds(1,3)+gmet(1,1)*(-6*gmet(2,3)&
&       *dgmetds(1,2)-3*gmet(2,2)*dgmetds(1,3))+gmet(1,3)*(-3*gmet(2,2)&
&       *dgmetds(1,1)+30*gmet(1,2)*dgmetds(1,2)-3*gmet(1,1)*dgmetds(2,2))&
&       -6*gmet(1,2)*(gmet(2,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(2,3))
       cm(10,1,2,3)=15*gmet(1,3)**2*dgmetds(1,3)-3*gmet(1,1)*gmet(3,3)&
&       *dgmetds(1,3)-3*gmet(1,3)*(gmet(3,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(3,3))
       cm(1,2,2,3)=9*gmet(1,2)**2*dgmetds(1,1)+18*gmet(1,1)*gmet(1,2)&
&       *dgmetds(1,2)+gmet(1,1)*(-6*gmet(2,2)*dgmetds(1,1)-3*gmet(1,1)&
&       *dgmetds(2,2))
       cm(2,2,2,3)=12*gmet(2,2)**2*dgmetds(1,1)+6*gmet(1,2)**2*dgmetds(2,2)&
&       +gmet(2,2)*(12*gmet(1,2)*dgmetds(1,2)+24*gmet(1,1)*dgmetds(2,2))
       cm(3,2,2,3)=15*gmet(2,3)**2*dgmetds(1,1)-12*gmet(1,2)*gmet(3,3)&
&       *dgmetds(1,2)-6*gmet(1,3)**2*dgmetds(2,2)-3*gmet(1,1)*gmet(3,3)&
&       *dgmetds(2,2)+18*gmet(1,2)*gmet(1,3)*dgmetds(2,3)+gmet(2,3)*(18*gmet(1,3)&
&       *dgmetds(1,2)+18*gmet(1,2)*dgmetds(1,3)+30*gmet(1,1)*dgmetds(2,3))&
&       -6*gmet(1,2)**2*dgmetds(3,3)+gmet(2,2)*(-3*gmet(3,3)*dgmetds(1,1)&
&       -12*gmet(1,3)*dgmetds(1,3)-3*gmet(1,1)*dgmetds(3,3))
       cm(4,2,2,3)=24*gmet(1,1)*gmet(2,3)*dgmetds(2,2)+gmet(1,2)*(12*gmet(2,3)&
&       *dgmetds(1,2)+6*gmet(1,3)*dgmetds(2,2))+6*gmet(1,2)**2*dgmetds(2,3)&
&       +gmet(2,2)*(24*gmet(2,3)*dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,2)&
&       +6*gmet(1,2)*dgmetds(1,3)+24*gmet(1,1)*dgmetds(2,3))
       cm(5,2,2,3)=3*gmet(1,2)**2*dgmetds(1,3)+gmet(1,1)*(24*gmet(2,3)&
&       *dgmetds(1,2)-9*gmet(2,2)*dgmetds(1,3))+gmet(1,3)*(-9*gmet(2,2)&
&       *dgmetds(1,1)+6*gmet(1,2)*dgmetds(1,2)-9*gmet(1,1)*dgmetds(2,2))&
&       +24*gmet(1,2)*(gmet(2,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(2,3))
       cm(6,2,2,3)=9*gmet(1,2)**2*dgmetds(1,2)+15*gmet(1,1)*gmet(2,2)&
&       *dgmetds(1,2)+15*gmet(1,2)*(gmet(2,2)*dgmetds(1,1)+gmet(1,1)&
&       *dgmetds(2,2))
       cm(7,2,2,3)=gmet(2,2)*(6*gmet(2,2)*dgmetds(1,2)+12*gmet(1,2)*dgmetds(2,2))
       cm(8,2,2,3)=3*gmet(2,3)**2*dgmetds(1,2)-9*gmet(1,2)*gmet(3,3)&
&       *dgmetds(2,2)+gmet(2,3)*(24*gmet(2,2)*dgmetds(1,3)+24*gmet(1,3)&
&       *dgmetds(2,2)+6*gmet(1,2)*dgmetds(2,3))+gmet(2,2)*(-9*gmet(3,3)&
&       *dgmetds(1,2)+24*gmet(1,3)*dgmetds(2,3)-9*gmet(1,2)*dgmetds(3,3))
       cm(9,2,2,3)=12*gmet(2,2)**2*dgmetds(1,3)+6*gmet(1,2)*gmet(2,3)&
&       *dgmetds(2,2)+gmet(2,2)*(6*gmet(2,3)*dgmetds(1,2)+24*gmet(1,3)&
&       *dgmetds(2,2)+6*gmet(1,2)*dgmetds(2,3))
       cm(10,2,2,3)=15*gmet(2,3)**2*dgmetds(1,3)+gmet(3,3)*(-3*gmet(1,3)&
&       *dgmetds(2,2)-6*gmet(1,2)*dgmetds(2,3))+gmet(2,3)*(-6*gmet(3,3)&
&       *dgmetds(1,2)+30*gmet(1,3)*dgmetds(2,3)-6*gmet(1,2)*dgmetds(3,3))&
&       -3*gmet(2,2)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(1,3,2,3)=9*gmet(1,3)**2*dgmetds(1,1)+18*gmet(1,1)*gmet(1,3)&
&       *dgmetds(1,3)+gmet(1,1)*(-6*gmet(3,3)*dgmetds(1,1)-3*gmet(1,1)&
&       *dgmetds(3,3))
       cm(2,3,2,3)=15*gmet(2,3)**2*dgmetds(1,1)-12*gmet(1,2)*gmet(3,3)&
&       *dgmetds(1,2)-6*gmet(1,3)**2*dgmetds(2,2)-3*gmet(1,1)*gmet(3,3)&
&       *dgmetds(2,2)+18*gmet(1,2)*gmet(1,3)*dgmetds(2,3)+gmet(2,3)*(18*gmet(1,3)&
&       *dgmetds(1,2)+18*gmet(1,2)*dgmetds(1,3)+30*gmet(1,1)*dgmetds(2,3))&
&       -6*gmet(1,2)**2*dgmetds(3,3)+gmet(2,2)*(-3*gmet(3,3)*dgmetds(1,1)&
&       -12*gmet(1,3)*dgmetds(1,3)-3*gmet(1,1)*dgmetds(3,3))
       cm(3,3,2,3)=12*gmet(3,3)**2*dgmetds(1,1)+6*gmet(1,3)**2*dgmetds(3,3)&
&       +gmet(3,3)*(12*gmet(1,3)*dgmetds(1,3)+24*gmet(1,1)*dgmetds(3,3))
       cm(4,3,2,3)=6*gmet(1,3)**2*dgmetds(2,3)+gmet(3,3)*(6*gmet(1,2)&
&       *dgmetds(1,3)+24*gmet(1,1)*dgmetds(2,3))+gmet(2,3)*(24*gmet(3,3)&
&       *dgmetds(1,1)+12*gmet(1,3)*dgmetds(1,3)+24*gmet(1,1)*dgmetds(3,3))&
&       +6*gmet(1,3)*(gmet(3,3)*dgmetds(1,2)+gmet(1,2)*dgmetds(3,3))
       cm(5,3,2,3)=9*gmet(1,3)**2*dgmetds(1,3)+15*gmet(1,1)*gmet(3,3)&
&       *dgmetds(1,3)+15*gmet(1,3)*(gmet(3,3)*dgmetds(1,1)+gmet(1,1)&
&       *dgmetds(3,3))
       cm(6,3,2,3)=3*gmet(1,3)**2*dgmetds(1,2)+gmet(1,1)*(-9*gmet(3,3)&
&       *dgmetds(1,2)+24*gmet(2,3)*dgmetds(1,3))+gmet(1,3)*(24*gmet(2,3)&
&       *dgmetds(1,1)+6*gmet(1,2)*dgmetds(1,3)+24*gmet(1,1)*dgmetds(2,3))&
&       -9*gmet(1,2)*(gmet(3,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(3,3))
       cm(7,3,2,3)=15*gmet(2,3)**2*dgmetds(1,2)-3*gmet(1,2)*gmet(3,3)&
&       *dgmetds(2,2)+gmet(2,3)*(-6*gmet(2,2)*dgmetds(1,3)-6*gmet(1,3)&
&       *dgmetds(2,2)+30*gmet(1,2)*dgmetds(2,3))+gmet(2,2)*(-3*gmet(3,3)&
&       *dgmetds(1,2)-6*gmet(1,3)*dgmetds(2,3)-3*gmet(1,2)*dgmetds(3,3))
       cm(8,3,2,3)=12*gmet(3,3)**2*dgmetds(1,2)+6*gmet(1,3)*gmet(2,3)&
&       *dgmetds(3,3)+gmet(3,3)*(6*gmet(2,3)*dgmetds(1,3)+6*gmet(1,3)&
&       *dgmetds(2,3)+24*gmet(1,2)*dgmetds(3,3))
       cm(9,3,2,3)=3*gmet(2,3)**2*dgmetds(1,3)+gmet(3,3)*(-9*gmet(1,3)&
&       *dgmetds(2,2)+24*gmet(1,2)*dgmetds(2,3))+gmet(2,3)*(24*gmet(3,3)&
&       *dgmetds(1,2)+6*gmet(1,3)*dgmetds(2,3)+24*gmet(1,2)*dgmetds(3,3))&
&       -9*gmet(2,2)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(10,3,2,3)=gmet(3,3)*(6*gmet(3,3)*dgmetds(1,3)+12*gmet(1,3)&
&       *dgmetds(3,3))
       cm(1,4,2,3)=18*gmet(1,2)*(gmet(1,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,3))&
&       +gmet(1,1)*(-12*gmet(2,3)*dgmetds(1,1)+18*gmet(1,3)*dgmetds(1,2)&
&       -6*gmet(1,1)*dgmetds(2,3))
       cm(2,4,2,3)=24*gmet(1,1)*gmet(2,3)*dgmetds(2,2)+gmet(1,2)*(12*gmet(2,3)&
&       *dgmetds(1,2)+6*gmet(1,3)*dgmetds(2,2))+6*gmet(1,2)**2*dgmetds(2,3)&
&       +gmet(2,2)*(24*gmet(2,3)*dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,2)&
&       +6*gmet(1,2)*dgmetds(1,3)+24*gmet(1,1)*dgmetds(2,3))
       cm(3,4,2,3)=6*gmet(1,3)**2*dgmetds(2,3)+gmet(3,3)*(6*gmet(1,2)&
&       *dgmetds(1,3)+24*gmet(1,1)*dgmetds(2,3))+gmet(2,3)*(24*gmet(3,3)&
&       *dgmetds(1,1)+12*gmet(1,3)*dgmetds(1,3)+24*gmet(1,1)*dgmetds(3,3))&
&       +6*gmet(1,3)*(gmet(3,3)*dgmetds(1,2)+gmet(1,2)*dgmetds(3,3))
       cm(4,4,2,3)=18*gmet(2,3)**2*dgmetds(1,1)+36*gmet(1,2)*gmet(3,3)&
&       *dgmetds(1,2)+18*gmet(1,3)**2*dgmetds(2,2)+30*gmet(1,1)*gmet(3,3)&
&       *dgmetds(2,2)-12*gmet(1,2)*gmet(1,3)*dgmetds(2,3)+gmet(2,3)*(-12*gmet(1,3)&
&       *dgmetds(1,2)-12*gmet(1,2)*dgmetds(1,3)+36*gmet(1,1)*dgmetds(2,3))&
&       +18*gmet(1,2)**2*dgmetds(3,3)+gmet(2,2)*(30*gmet(3,3)*dgmetds(1,1)&
&       +36*gmet(1,3)*dgmetds(1,3)+30*gmet(1,1)*dgmetds(3,3))
       cm(5,4,2,3)=6*gmet(1,3)**2*dgmetds(1,2)+gmet(1,1)*(24*gmet(3,3)&
&       *dgmetds(1,2)+6*gmet(2,3)*dgmetds(1,3))+gmet(1,3)*(6*gmet(2,3)&
&       *dgmetds(1,1)+12*gmet(1,2)*dgmetds(1,3)+6*gmet(1,1)*dgmetds(2,3))&
&       +24*gmet(1,2)*(gmet(3,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(3,3))
       cm(6,4,2,3)=6*gmet(1,2)**2*dgmetds(1,3)+gmet(1,1)*(6*gmet(2,3)&
&       *dgmetds(1,2)+24*gmet(2,2)*dgmetds(1,3))+gmet(1,3)*(24*gmet(2,2)&
&       *dgmetds(1,1)+12*gmet(1,2)*dgmetds(1,2)+24*gmet(1,1)*dgmetds(2,2))&
&       +6*gmet(1,2)*(gmet(2,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(2,3))
       cm(7,4,2,3)=-6*gmet(2,2)**2*dgmetds(1,3)+18*gmet(1,2)*gmet(2,3)&
&       *dgmetds(2,2)+gmet(2,2)*(18*gmet(2,3)*dgmetds(1,2)-12*gmet(1,3)&
&       *dgmetds(2,2)+18*gmet(1,2)*dgmetds(2,3))
       cm(8,4,2,3)=6*gmet(2,3)**2*dgmetds(1,3)+gmet(3,3)*(24*gmet(1,3)&
&       *dgmetds(2,2)+6*gmet(1,2)*dgmetds(2,3))+gmet(2,3)*(6*gmet(3,3)&
&       *dgmetds(1,2)+12*gmet(1,3)*dgmetds(2,3)+6*gmet(1,2)*dgmetds(3,3))&
&       +24*gmet(2,2)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(9,4,2,3)=6*gmet(2,3)**2*dgmetds(1,2)+24*gmet(1,2)*gmet(3,3)&
&       *dgmetds(2,2)+gmet(2,3)*(6*gmet(2,2)*dgmetds(1,3)+6*gmet(1,3)&
&       *dgmetds(2,2)+12*gmet(1,2)*dgmetds(2,3))+gmet(2,2)*(24*gmet(3,3)&
&       *dgmetds(1,2)+6*gmet(1,3)*dgmetds(2,3)+24*gmet(1,2)*dgmetds(3,3))
       cm(10,4,2,3)=-6*gmet(3,3)**2*dgmetds(1,2)+18*gmet(1,3)*gmet(2,3)&
&       *dgmetds(3,3)+gmet(3,3)*(18*gmet(2,3)*dgmetds(1,3)+18*gmet(1,3)&
&       *dgmetds(2,3)-12*gmet(1,2)*dgmetds(3,3))
       cm(1,5,2,3)=gmet(1,1)*(12*gmet(1,3)*dgmetds(1,1)+6*gmet(1,1)*dgmetds(1,3))
       cm(2,5,2,3)=3*gmet(1,2)**2*dgmetds(1,3)+gmet(1,1)*(24*gmet(2,3)&
&       *dgmetds(1,2)-9*gmet(2,2)*dgmetds(1,3))+gmet(1,3)*(-9*gmet(2,2)&
&       *dgmetds(1,1)+6*gmet(1,2)*dgmetds(1,2)-9*gmet(1,1)*dgmetds(2,2))&
&       +24*gmet(1,2)*(gmet(2,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(2,3))
       cm(3,5,2,3)=9*gmet(1,3)**2*dgmetds(1,3)+15*gmet(1,1)*gmet(3,3)&
&       *dgmetds(1,3)+15*gmet(1,3)*(gmet(3,3)*dgmetds(1,1)+gmet(1,1)&
&       *dgmetds(3,3))
       cm(4,5,2,3)=6*gmet(1,3)**2*dgmetds(1,2)+gmet(1,1)*(24*gmet(3,3)&
&       *dgmetds(1,2)+6*gmet(2,3)*dgmetds(1,3))+gmet(1,3)*(6*gmet(2,3)&
&       *dgmetds(1,1)+12*gmet(1,2)*dgmetds(1,3)+6*gmet(1,1)*dgmetds(2,3))&
&       +24*gmet(1,2)*(gmet(3,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(3,3))
       cm(5,5,2,3)=6*gmet(1,3)**2*dgmetds(1,1)+12*gmet(1,1)*gmet(1,3)&
&       *dgmetds(1,3)+gmet(1,1)*(24*gmet(3,3)*dgmetds(1,1)+12*gmet(1,1)&
&       *dgmetds(3,3))
       cm(6,5,2,3)=6*gmet(1,2)*(gmet(1,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,3))&
&       +gmet(1,1)*(24*gmet(2,3)*dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,2)&
&       +12*gmet(1,1)*dgmetds(2,3))
       cm(7,5,2,3)=-3*gmet(1,1)*gmet(2,3)*dgmetds(2,2)+gmet(1,2)*(30*gmet(2,3)&
&       *dgmetds(1,2)-6*gmet(1,3)*dgmetds(2,2))+15*gmet(1,2)**2*dgmetds(2,3)&
&       +gmet(2,2)*(-3*gmet(2,3)*dgmetds(1,1)-6*gmet(1,3)*dgmetds(1,2)&
&       -6*gmet(1,2)*dgmetds(1,3)-3*gmet(1,1)*dgmetds(2,3))
       cm(8,5,2,3)=3*gmet(1,3)**2*dgmetds(2,3)+gmet(3,3)*(24*gmet(1,2)&
&       *dgmetds(1,3)-9*gmet(1,1)*dgmetds(2,3))+gmet(2,3)*(-9*gmet(3,3)&
&       *dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,3)-9*gmet(1,1)*dgmetds(3,3))&
&       +24*gmet(1,3)*(gmet(3,3)*dgmetds(1,2)+gmet(1,2)*dgmetds(3,3))
       cm(9,5,2,3)=-6*gmet(2,3)**2*dgmetds(1,1)+30*gmet(1,2)*gmet(3,3)&
&       *dgmetds(1,2)-6*gmet(1,3)**2*dgmetds(2,2)-3*gmet(1,1)*gmet(3,3)&
&       *dgmetds(2,2)+18*gmet(1,2)*gmet(1,3)*dgmetds(2,3)+gmet(2,3)*(18*gmet(1,3)&
&       *dgmetds(1,2)+18*gmet(1,2)*dgmetds(1,3)-12*gmet(1,1)*dgmetds(2,3))&
&       +15*gmet(1,2)**2*dgmetds(3,3)+gmet(2,2)*(-3*gmet(3,3)*dgmetds(1,1)&
&       -12*gmet(1,3)*dgmetds(1,3)-3*gmet(1,1)*dgmetds(3,3))
       cm(10,5,2,3)=-3*gmet(3,3)**2*dgmetds(1,1)+9*gmet(1,3)**2*dgmetds(3,3)&
&       +gmet(3,3)*(18*gmet(1,3)*dgmetds(1,3)-6*gmet(1,1)*dgmetds(3,3))
       cm(1,6,2,3)=gmet(1,1)*(12*gmet(1,2)*dgmetds(1,1)+6*gmet(1,1)*dgmetds(1,2))
       cm(2,6,2,3)=9*gmet(1,2)**2*dgmetds(1,2)+15*gmet(1,1)*gmet(2,2)&
&       *dgmetds(1,2)+15*gmet(1,2)*(gmet(2,2)*dgmetds(1,1)+gmet(1,1)&
&       *dgmetds(2,2))
       cm(3,6,2,3)=3*gmet(1,3)**2*dgmetds(1,2)+gmet(1,1)*(-9*gmet(3,3)&
&       *dgmetds(1,2)+24*gmet(2,3)*dgmetds(1,3))+gmet(1,3)*(24*gmet(2,3)&
&       *dgmetds(1,1)+6*gmet(1,2)*dgmetds(1,3)+24*gmet(1,1)*dgmetds(2,3))&
&       -9*gmet(1,2)*(gmet(3,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(3,3))
       cm(4,6,2,3)=6*gmet(1,2)**2*dgmetds(1,3)+gmet(1,1)*(6*gmet(2,3)&
&       *dgmetds(1,2)+24*gmet(2,2)*dgmetds(1,3))+gmet(1,3)*(24*gmet(2,2)&
&       *dgmetds(1,1)+12*gmet(1,2)*dgmetds(1,2)+24*gmet(1,1)*dgmetds(2,2))&
&       +6*gmet(1,2)*(gmet(2,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(2,3))
       cm(5,6,2,3)=6*gmet(1,2)*(gmet(1,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,3))&
&       +gmet(1,1)*(24*gmet(2,3)*dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,2)&
&       +12*gmet(1,1)*dgmetds(2,3))
       cm(6,6,2,3)=6*gmet(1,2)**2*dgmetds(1,1)+12*gmet(1,1)*gmet(1,2)&
&       *dgmetds(1,2)+gmet(1,1)*(24*gmet(2,2)*dgmetds(1,1)+12*gmet(1,1)&
&       *dgmetds(2,2))
       cm(7,6,2,3)=-3*gmet(2,2)**2*dgmetds(1,1)+9*gmet(1,2)**2*dgmetds(2,2)&
&       +gmet(2,2)*(18*gmet(1,2)*dgmetds(1,2)-6*gmet(1,1)*dgmetds(2,2))
       cm(8,6,2,3)=-6*gmet(2,3)**2*dgmetds(1,1)-12*gmet(1,2)*gmet(3,3)&
&       *dgmetds(1,2)+15*gmet(1,3)**2*dgmetds(2,2)-3*gmet(1,1)*gmet(3,3)&
&       *dgmetds(2,2)+18*gmet(1,2)*gmet(1,3)*dgmetds(2,3)+gmet(2,3)*(18*gmet(1,3)&
&       *dgmetds(1,2)+18*gmet(1,2)*dgmetds(1,3)-12*gmet(1,1)*dgmetds(2,3))&
&       -6*gmet(1,2)**2*dgmetds(3,3)+gmet(2,2)*(-3*gmet(3,3)*dgmetds(1,1)&
&       +30*gmet(1,3)*dgmetds(1,3)-3*gmet(1,1)*dgmetds(3,3))
       cm(9,6,2,3)=-9*gmet(1,1)*gmet(2,3)*dgmetds(2,2)+gmet(1,2)*(6*gmet(2,3)&
&       *dgmetds(1,2)+24*gmet(1,3)*dgmetds(2,2))+3*gmet(1,2)**2*dgmetds(2,3)&
&       +gmet(2,2)*(-9*gmet(2,3)*dgmetds(1,1)+24*gmet(1,3)*dgmetds(1,2)&
&       +24*gmet(1,2)*dgmetds(1,3)-9*gmet(1,1)*dgmetds(2,3))
       cm(10,6,2,3)=15*gmet(1,3)**2*dgmetds(2,3)+gmet(3,3)*(-6*gmet(1,2)&
&       *dgmetds(1,3)-3*gmet(1,1)*dgmetds(2,3))+gmet(2,3)*(-3*gmet(3,3)&
&       *dgmetds(1,1)+30*gmet(1,3)*dgmetds(1,3)-3*gmet(1,1)*dgmetds(3,3))&
&       -6*gmet(1,3)*(gmet(3,3)*dgmetds(1,2)+gmet(1,2)*dgmetds(3,3))
       cm(1,7,2,3)=15*gmet(1,2)**2*dgmetds(1,2)-3*gmet(1,1)*gmet(2,2)&
&       *dgmetds(1,2)-3*gmet(1,2)*(gmet(2,2)*dgmetds(1,1)+gmet(1,1)*dgmetds(2,2))
       cm(2,7,2,3)=gmet(2,2)*(6*gmet(2,2)*dgmetds(1,2)+12*gmet(1,2)*dgmetds(2,2))
       cm(3,7,2,3)=15*gmet(2,3)**2*dgmetds(1,2)-3*gmet(1,2)*gmet(3,3)&
&       *dgmetds(2,2)+gmet(2,3)*(-6*gmet(2,2)*dgmetds(1,3)-6*gmet(1,3)&
&       *dgmetds(2,2)+30*gmet(1,2)*dgmetds(2,3))+gmet(2,2)*(-3*gmet(3,3)&
&       *dgmetds(1,2)-6*gmet(1,3)*dgmetds(2,3)-3*gmet(1,2)*dgmetds(3,3))
       cm(4,7,2,3)=-6*gmet(2,2)**2*dgmetds(1,3)+18*gmet(1,2)*gmet(2,3)&
&       *dgmetds(2,2)+gmet(2,2)*(18*gmet(2,3)*dgmetds(1,2)-12*gmet(1,3)&
&       *dgmetds(2,2)+18*gmet(1,2)*dgmetds(2,3))
       cm(5,7,2,3)=-3*gmet(1,1)*gmet(2,3)*dgmetds(2,2)+gmet(1,2)*(30*gmet(2,3)&
&       *dgmetds(1,2)-6*gmet(1,3)*dgmetds(2,2))+15*gmet(1,2)**2*dgmetds(2,3)&
&       +gmet(2,2)*(-3*gmet(2,3)*dgmetds(1,1)-6*gmet(1,3)*dgmetds(1,2)&
&       -6*gmet(1,2)*dgmetds(1,3)-3*gmet(1,1)*dgmetds(2,3))
       cm(6,7,2,3)=-3*gmet(2,2)**2*dgmetds(1,1)+9*gmet(1,2)**2*dgmetds(2,2)&
&       +gmet(2,2)*(18*gmet(1,2)*dgmetds(1,2)-6*gmet(1,1)*dgmetds(2,2))
       cm(7,7,2,3)=6*gmet(2,2)**2*dgmetds(2,2)
       cm(8,7,2,3)=9*gmet(2,3)**2*dgmetds(2,2)+18*gmet(2,2)*gmet(2,3)&
&       *dgmetds(2,3)+gmet(2,2)*(-6*gmet(3,3)*dgmetds(2,2)-3*gmet(2,2)&
&       *dgmetds(3,3))
       cm(9,7,2,3)=gmet(2,2)*(12*gmet(2,3)*dgmetds(2,2)+6*gmet(2,2)*dgmetds(2,3))
       cm(10,7,2,3)=15*gmet(2,3)**2*dgmetds(2,3)-3*gmet(2,2)*gmet(3,3)&
&       *dgmetds(2,3)-3*gmet(2,3)*(gmet(3,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(3,3))
       cm(1,8,2,3)=15*gmet(1,3)**2*dgmetds(1,2)+gmet(1,1)*(-3*gmet(3,3)&
&       *dgmetds(1,2)-6*gmet(2,3)*dgmetds(1,3))+gmet(1,3)*(-6*gmet(2,3)&
&       *dgmetds(1,1)+30*gmet(1,2)*dgmetds(1,3)-6*gmet(1,1)*dgmetds(2,3))&
&       -3*gmet(1,2)*(gmet(3,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(3,3))
       cm(2,8,2,3)=3*gmet(2,3)**2*dgmetds(1,2)-9*gmet(1,2)*gmet(3,3)&
&       *dgmetds(2,2)+gmet(2,3)*(24*gmet(2,2)*dgmetds(1,3)+24*gmet(1,3)&
&       *dgmetds(2,2)+6*gmet(1,2)*dgmetds(2,3))+gmet(2,2)*(-9*gmet(3,3)&
&       *dgmetds(1,2)+24*gmet(1,3)*dgmetds(2,3)-9*gmet(1,2)*dgmetds(3,3))
       cm(3,8,2,3)=12*gmet(3,3)**2*dgmetds(1,2)+6*gmet(1,3)*gmet(2,3)&
&       *dgmetds(3,3)+gmet(3,3)*(6*gmet(2,3)*dgmetds(1,3)+6*gmet(1,3)&
&       *dgmetds(2,3)+24*gmet(1,2)*dgmetds(3,3))
       cm(4,8,2,3)=6*gmet(2,3)**2*dgmetds(1,3)+gmet(3,3)*(24*gmet(1,3)&
&       *dgmetds(2,2)+6*gmet(1,2)*dgmetds(2,3))+gmet(2,3)*(6*gmet(3,3)&
&       *dgmetds(1,2)+12*gmet(1,3)*dgmetds(2,3)+6*gmet(1,2)*dgmetds(3,3))&
&       +24*gmet(2,2)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(5,8,2,3)=3*gmet(1,3)**2*dgmetds(2,3)+gmet(3,3)*(24*gmet(1,2)&
&       *dgmetds(1,3)-9*gmet(1,1)*dgmetds(2,3))+gmet(2,3)*(-9*gmet(3,3)&
&       *dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,3)-9*gmet(1,1)*dgmetds(3,3))&
&       +24*gmet(1,3)*(gmet(3,3)*dgmetds(1,2)+gmet(1,2)*dgmetds(3,3))
       cm(6,8,2,3)=-6*gmet(2,3)**2*dgmetds(1,1)-12*gmet(1,2)*gmet(3,3)&
&       *dgmetds(1,2)+15*gmet(1,3)**2*dgmetds(2,2)-3*gmet(1,1)*gmet(3,3)&
&       *dgmetds(2,2)+18*gmet(1,2)*gmet(1,3)*dgmetds(2,3)+gmet(2,3)*(18*gmet(1,3)&
&       *dgmetds(1,2)+18*gmet(1,2)*dgmetds(1,3)-12*gmet(1,1)*dgmetds(2,3))&
&       -6*gmet(1,2)**2*dgmetds(3,3)+gmet(2,2)*(-3*gmet(3,3)*dgmetds(1,1)&
&       +30*gmet(1,3)*dgmetds(1,3)-3*gmet(1,1)*dgmetds(3,3))
       cm(7,8,2,3)=9*gmet(2,3)**2*dgmetds(2,2)+18*gmet(2,2)*gmet(2,3)&
&       *dgmetds(2,3)+gmet(2,2)*(-6*gmet(3,3)*dgmetds(2,2)-3*gmet(2,2)&
&       *dgmetds(3,3))
       cm(8,8,2,3)=12*gmet(3,3)**2*dgmetds(2,2)+6*gmet(2,3)**2*dgmetds(3,3)&
&       +gmet(3,3)*(12*gmet(2,3)*dgmetds(2,3)+24*gmet(2,2)*dgmetds(3,3))
       cm(9,8,2,3)=9*gmet(2,3)**2*dgmetds(2,3)+15*gmet(2,2)*gmet(3,3)&
&       *dgmetds(2,3)+15*gmet(2,3)*(gmet(3,3)*dgmetds(2,2)+gmet(2,2)&
&       *dgmetds(3,3))
       cm(10,8,2,3)=gmet(3,3)*(6*gmet(3,3)*dgmetds(2,3)+12*gmet(2,3)&
&       *dgmetds(3,3))
       cm(1,9,2,3)=15*gmet(1,2)**2*dgmetds(1,3)+gmet(1,1)*(-6*gmet(2,3)&
&       *dgmetds(1,2)-3*gmet(2,2)*dgmetds(1,3))+gmet(1,3)*(-3*gmet(2,2)&
&       *dgmetds(1,1)+30*gmet(1,2)*dgmetds(1,2)-3*gmet(1,1)*dgmetds(2,2))&
&       -6*gmet(1,2)*(gmet(2,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(2,3))
       cm(2,9,2,3)=12*gmet(2,2)**2*dgmetds(1,3)+6*gmet(1,2)*gmet(2,3)&
&       *dgmetds(2,2)+gmet(2,2)*(6*gmet(2,3)*dgmetds(1,2)+24*gmet(1,3)&
&       *dgmetds(2,2)+6*gmet(1,2)*dgmetds(2,3))
       cm(3,9,2,3)=3*gmet(2,3)**2*dgmetds(1,3)+gmet(3,3)*(-9*gmet(1,3)&
&       *dgmetds(2,2)+24*gmet(1,2)*dgmetds(2,3))+gmet(2,3)*(24*gmet(3,3)&
&       *dgmetds(1,2)+6*gmet(1,3)*dgmetds(2,3)+24*gmet(1,2)*dgmetds(3,3))&
&       -9*gmet(2,2)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(4,9,2,3)=6*gmet(2,3)**2*dgmetds(1,2)+24*gmet(1,2)*gmet(3,3)&
&       *dgmetds(2,2)+gmet(2,3)*(6*gmet(2,2)*dgmetds(1,3)+6*gmet(1,3)&
&       *dgmetds(2,2)+12*gmet(1,2)*dgmetds(2,3))+gmet(2,2)*(24*gmet(3,3)&
&       *dgmetds(1,2)+6*gmet(1,3)*dgmetds(2,3)+24*gmet(1,2)*dgmetds(3,3))
       cm(5,9,2,3)=-6*gmet(2,3)**2*dgmetds(1,1)+30*gmet(1,2)*gmet(3,3)&
&       *dgmetds(1,2)-6*gmet(1,3)**2*dgmetds(2,2)-3*gmet(1,1)*gmet(3,3)&
&       *dgmetds(2,2)+18*gmet(1,2)*gmet(1,3)*dgmetds(2,3)+gmet(2,3)*(18*gmet(1,3)&
&       *dgmetds(1,2)+18*gmet(1,2)*dgmetds(1,3)-12*gmet(1,1)*dgmetds(2,3))&
&       +15*gmet(1,2)**2*dgmetds(3,3)+gmet(2,2)*(-3*gmet(3,3)*dgmetds(1,1)&
&       -12*gmet(1,3)*dgmetds(1,3)-3*gmet(1,1)*dgmetds(3,3))
       cm(6,9,2,3)=-9*gmet(1,1)*gmet(2,3)*dgmetds(2,2)+gmet(1,2)*(6*gmet(2,3)&
&       *dgmetds(1,2)+24*gmet(1,3)*dgmetds(2,2))+3*gmet(1,2)**2*dgmetds(2,3)&
&       +gmet(2,2)*(-9*gmet(2,3)*dgmetds(1,1)+24*gmet(1,3)*dgmetds(1,2)&
&       +24*gmet(1,2)*dgmetds(1,3)-9*gmet(1,1)*dgmetds(2,3))
       cm(7,9,2,3)=gmet(2,2)*(12*gmet(2,3)*dgmetds(2,2)+6*gmet(2,2)*dgmetds(2,3))
       cm(8,9,2,3)=9*gmet(2,3)**2*dgmetds(2,3)+15*gmet(2,2)*gmet(3,3)&
&       *dgmetds(2,3)+15*gmet(2,3)*(gmet(3,3)*dgmetds(2,2)+gmet(2,2)&
&       *dgmetds(3,3))
       cm(9,9,2,3)=6*gmet(2,3)**2*dgmetds(2,2)+12*gmet(2,2)*gmet(2,3)&
&       *dgmetds(2,3)+gmet(2,2)*(24*gmet(3,3)*dgmetds(2,2)+12*gmet(2,2)&
&       *dgmetds(3,3))
       cm(10,9,2,3)=-3*gmet(3,3)**2*dgmetds(2,2)+9*gmet(2,3)**2*dgmetds(3,3)&
&       +gmet(3,3)*(18*gmet(2,3)*dgmetds(2,3)-6*gmet(2,2)*dgmetds(3,3))
       cm(1,10,2,3)=15*gmet(1,3)**2*dgmetds(1,3)-3*gmet(1,1)*gmet(3,3)&
&       *dgmetds(1,3)-3*gmet(1,3)*(gmet(3,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(3,3))
       cm(2,10,2,3)=15*gmet(2,3)**2*dgmetds(1,3)+gmet(3,3)*(-3*gmet(1,3)&
&       *dgmetds(2,2)-6*gmet(1,2)*dgmetds(2,3))+gmet(2,3)*(-6*gmet(3,3)&
&       *dgmetds(1,2)+30*gmet(1,3)*dgmetds(2,3)-6*gmet(1,2)*dgmetds(3,3))&
&       -3*gmet(2,2)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(3,10,2,3)=gmet(3,3)*(6*gmet(3,3)*dgmetds(1,3)+12*gmet(1,3)&
&       *dgmetds(3,3))
       cm(4,10,2,3)=-6*gmet(3,3)**2*dgmetds(1,2)+18*gmet(1,3)*gmet(2,3)&
&       *dgmetds(3,3)+gmet(3,3)*(18*gmet(2,3)*dgmetds(1,3)+18*gmet(1,3)&
&       *dgmetds(2,3)-12*gmet(1,2)*dgmetds(3,3))
       cm(5,10,2,3)=-3*gmet(3,3)**2*dgmetds(1,1)+9*gmet(1,3)**2*dgmetds(3,3)&
&       +gmet(3,3)*(18*gmet(1,3)*dgmetds(1,3)-6*gmet(1,1)*dgmetds(3,3))
       cm(6,10,2,3)=15*gmet(1,3)**2*dgmetds(2,3)+gmet(3,3)*(-6*gmet(1,2)&
&       *dgmetds(1,3)-3*gmet(1,1)*dgmetds(2,3))+gmet(2,3)*(-3*gmet(3,3)&
&       *dgmetds(1,1)+30*gmet(1,3)*dgmetds(1,3)-3*gmet(1,1)*dgmetds(3,3))&
&       -6*gmet(1,3)*(gmet(3,3)*dgmetds(1,2)+gmet(1,2)*dgmetds(3,3))
       cm(7,10,2,3)=15*gmet(2,3)**2*dgmetds(2,3)-3*gmet(2,2)*gmet(3,3)&
&       *dgmetds(2,3)-3*gmet(2,3)*(gmet(3,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(3,3))
       cm(8,10,2,3)=gmet(3,3)*(6*gmet(3,3)*dgmetds(2,3)+12*gmet(2,3)&
&       *dgmetds(3,3))
       cm(9,10,2,3)=-3*gmet(3,3)**2*dgmetds(2,2)+9*gmet(2,3)**2*dgmetds(3,3)&
&       +gmet(3,3)*(18*gmet(2,3)*dgmetds(2,3)-6*gmet(2,2)*dgmetds(3,3))
       cm(10,10,2,3)=6*gmet(3,3)**2*dgmetds(3,3)
     elseif(iterm==3)then
       cm(1,1,3,3)=gmet(1,1)**3*dgmetds(1,1)
       cm(2,1,3,3)=gmet(1,1)*(4.5d0*gmet(1,2)**2*dgmetds(1,1)+6*gmet(1,1)&
&       *gmet(1,2)*dgmetds(1,2)+gmet(1,1)*(-1.5d0*gmet(2,2)*dgmetds(1,1)&
&       +gmet(1,1)*dgmetds(2,2)))
       cm(3,1,3,3)=gmet(1,1)*(4.5d0*gmet(1,3)**2*dgmetds(1,1)+6*gmet(1,1)&
&       *gmet(1,3)*dgmetds(1,3)+gmet(1,1)*(-1.5d0*gmet(3,3)*dgmetds(1,1)&
&       +gmet(1,1)*dgmetds(3,3)))
       cm(4,1,3,3)=gmet(1,1)*(gmet(1,2)*(9*gmet(1,3)*dgmetds(1,1)+6*gmet(1,1)&
&       *dgmetds(1,3))+gmet(1,1)*(-3*gmet(2,3)*dgmetds(1,1)+6*gmet(1,3)&
&       *dgmetds(1,2)+2*gmet(1,1)*dgmetds(2,3)))
       cm(5,1,3,3)=gmet(1,1)**2*(3*gmet(1,3)*dgmetds(1,1)+2*gmet(1,1)&
&       *dgmetds(1,3))
       cm(6,1,3,3)=gmet(1,1)**2*(3*gmet(1,2)*dgmetds(1,1)+2*gmet(1,1)&
&       *dgmetds(1,2))
       cm(7,1,3,3)=2.5d0*gmet(1,2)**3*dgmetds(1,1)+9*gmet(1,1)*gmet(1,2)&
&       **2*dgmetds(1,2)-3*gmet(1,1)**2*gmet(2,2)*dgmetds(1,2)+gmet(1,1)&
&       *gmet(1,2)*(-1.5d0*gmet(2,2)*dgmetds(1,1)+3*gmet(1,1)*dgmetds(2,2))
       cm(8,1,3,3)=(gmet(1,1)*(216*gmet(1,3)**2*dgmetds(1,2)+gmet(1,1)&
&       *(-72*gmet(3,3)*dgmetds(1,2)-144*gmet(2,3)*dgmetds(1,3))+gmet(1,3)&
&       *(-72*gmet(2,3)*dgmetds(1,1)+144*gmet(1,1)*dgmetds(2,3)))+gmet(1,2)&
&       *(180*gmet(1,3)**2*dgmetds(1,1)+432*gmet(1,1)*gmet(1,3)*dgmetds(1,3)&
&       +gmet(1,1)*(-36*gmet(3,3)*dgmetds(1,1)+72*gmet(1,1)*dgmetds(3,3))))&
&       /24.d0
       cm(9,1,3,3)=(gmet(1,2)**2*(180*gmet(1,3)*dgmetds(1,1)+216*gmet(1,1)&
&       *dgmetds(1,3))+gmet(1,1)*(gmet(1,1)*(-144*gmet(2,3)*dgmetds(1,2)&
&       -72*gmet(2,2)*dgmetds(1,3))+gmet(1,3)*(-36*gmet(2,2)*dgmetds(1,1)&
&       +72*gmet(1,1)*dgmetds(2,2)))+gmet(1,1)*gmet(1,2)*(-72*gmet(2,3)&
&       *dgmetds(1,1)+432*gmet(1,3)*dgmetds(1,2)+144*gmet(1,1)*dgmetds(2,3)))&
&       /24.d0
       cm(10,1,3,3)=2.5d0*gmet(1,3)**3*dgmetds(1,1)+9*gmet(1,1)*gmet(1,3)&
&       **2*dgmetds(1,3)-3*gmet(1,1)**2*gmet(3,3)*dgmetds(1,3)+gmet(1,1)&
&       *gmet(1,3)*(-1.5d0*gmet(3,3)*dgmetds(1,1)+3*gmet(1,1)*dgmetds(3,3))
       cm(11,1,3,3)=5*gmet(1,2)**3*dgmetds(1,2)-3*gmet(1,1)*gmet(1,2)&
&       *gmet(2,2)*dgmetds(1,2)+4.5d0*gmet(1,1)*gmet(1,2)**2*dgmetds(2,2)&
&       -1.5d0*gmet(1,1)**2*gmet(2,2)*dgmetds(2,2)
       cm(12,1,3,3)=(2*(-36*gmet(1,1)*gmet(1,3)*gmet(2,3)+gmet(1,2)*(90*gmet(1,3)&
&       **2-18*gmet(1,1)*gmet(3,3)))*dgmetds(1,2)+2*(90*gmet(1,2)**2*gmet(1,3)&
&       -18*gmet(1,1)*gmet(1,3)*gmet(2,2)-36*gmet(1,1)*gmet(1,2)*gmet(2,3))&
&       *dgmetds(1,3)+gmet(1,1)*(54*gmet(1,3)**2-18*gmet(1,1)*gmet(3,3))&
&       *dgmetds(2,2)+4*gmet(1,1)*(54*gmet(1,2)*gmet(1,3)-18*gmet(1,1)&
&       *gmet(2,3))*dgmetds(2,3)+gmet(1,1)*(54*gmet(1,2)**2-18*gmet(1,1)&
&       *gmet(2,2))*dgmetds(3,3))/12.d0
       cm(13,1,3,3)=(180*gmet(1,2)**3*dgmetds(1,3)+gmet(1,1)*gmet(1,2)&
&       *(-216*gmet(2,3)*dgmetds(1,2)-108*gmet(2,2)*dgmetds(1,3)+324*gmet(1,3)&
&       *dgmetds(2,2))+gmet(1,2)**2*(540*gmet(1,3)*dgmetds(1,2)+324*gmet(1,1)&
&       *dgmetds(2,3))-108*gmet(1,1)*(gmet(1,3)*gmet(2,2)*dgmetds(1,2)&
&       +gmet(1,1)*(gmet(2,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(2,3))))&
&       /36.d0
       cm(14,1,3,3)=(180*gmet(1,3)**3*dgmetds(1,2)+gmet(1,3)**2*(540*gmet(1,2)&
&       *dgmetds(1,3)+324*gmet(1,1)*dgmetds(2,3))+gmet(1,1)*gmet(1,3)&
&       *(-108*gmet(3,3)*dgmetds(1,2)-216*gmet(2,3)*dgmetds(1,3)+324*gmet(1,2)&
&       *dgmetds(3,3))-108*gmet(1,1)*(gmet(1,2)*gmet(3,3)*dgmetds(1,3)&
&       +gmet(1,1)*(gmet(3,3)*dgmetds(2,3)+gmet(2,3)*dgmetds(3,3))))&
&       /36.d0
       cm(15,1,3,3)=5*gmet(1,3)**3*dgmetds(1,3)-3*gmet(1,1)*gmet(1,3)&
&       *gmet(3,3)*dgmetds(1,3)+4.5d0*gmet(1,1)*gmet(1,3)**2*dgmetds(3,3)&
&       -1.5d0*gmet(1,1)**2*gmet(3,3)*dgmetds(3,3)
       cm(16,1,3,3)=gmet(1,2)*(2.5d0*gmet(1,2)**2-1.5d0*gmet(1,1)*gmet(2,2))&
&       *dgmetds(2,2)
       cm(17,1,3,3)=(1080*gmet(1,2)**2*gmet(1,3)*dgmetds(2,3)-216*gmet(1,1)&
&       *gmet(1,3)*(gmet(2,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(2,3))+180*gmet(1,2)&
&       **3*dgmetds(3,3)+gmet(1,2)*(540*gmet(1,3)**2*dgmetds(2,2)+gmet(1,1)&
&       *(-108*gmet(3,3)*dgmetds(2,2)-432*gmet(2,3)*dgmetds(2,3)-108*gmet(2,2)&
&       *dgmetds(3,3))))/72.d0
       cm(18,1,3,3)=7.5d0*gmet(1,2)**2*gmet(1,3)*dgmetds(2,2)-1.5d0*gmet(1,1)&
&       *gmet(1,3)*gmet(2,2)*dgmetds(2,2)+5*gmet(1,2)**3*dgmetds(2,3)&
&       -3*gmet(1,1)*gmet(1,2)*(gmet(2,3)*dgmetds(2,2)+gmet(2,2)*dgmetds(2,3))
       cm(19,1,3,3)=(180*gmet(1,3)**3*dgmetds(2,2)+1080*gmet(1,2)*gmet(1,3)&
&       **2*dgmetds(2,3)-216*gmet(1,1)*gmet(1,2)*(gmet(3,3)*dgmetds(2,3)&
&       +gmet(2,3)*dgmetds(3,3))+gmet(1,3)*(540*gmet(1,2)**2*dgmetds(3,3)&
&       +gmet(1,1)*(-108*gmet(3,3)*dgmetds(2,2)-432*gmet(2,3)*dgmetds(2,3)&
&       -108*gmet(2,2)*dgmetds(3,3))))/72.d0
       cm(20,1,3,3)=5*gmet(1,3)**3*dgmetds(2,3)+7.5d0*gmet(1,2)*gmet(1,3)&
&       **2*dgmetds(3,3)-1.5d0*gmet(1,1)*gmet(1,2)*gmet(3,3)*dgmetds(3,3)&
&       -3*gmet(1,1)*gmet(1,3)*(gmet(3,3)*dgmetds(2,3)+gmet(2,3)*dgmetds(3,3))
       cm(21,1,3,3)=gmet(1,3)*(2.5d0*gmet(1,3)**2-1.5d0*gmet(1,1)*gmet(3,3))&
&       *dgmetds(3,3)
       cm(1,2,3,3)=gmet(1,1)*(4.5d0*gmet(1,2)**2-1.5d0*gmet(1,1)*gmet(2,2))&
&       *dgmetds(1,1)
       cm(2,2,3,3)=3*gmet(1,2)**3*dgmetds(1,2)+15*gmet(1,1)*gmet(1,2)&
&       *gmet(2,2)*dgmetds(1,2)+gmet(1,1)*gmet(2,2)*(6*gmet(2,2)*dgmetds(1,1)&
&       -1.5d0*gmet(1,1)*dgmetds(2,2))+gmet(1,2)**2*(3*gmet(2,2)*dgmetds(1,1)&
&       +4.5d0*gmet(1,1)*dgmetds(2,2))
       cm(3,2,3,3)=(6*(-12*gmet(1,3)**2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(30*gmet(2,3)&
&       **2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+12*(6*gmet(1,2)**2*gmet(1,3)&
&       -18*gmet(1,1)*gmet(1,3)*gmet(2,2)+48*gmet(1,1)*gmet(1,2)*gmet(2,3))&
&       *dgmetds(1,3)+2*gmet(1,1)*(54*gmet(1,2)**2-18*gmet(1,1)*gmet(2,2))&
&       *dgmetds(3,3))/24.d0
       cm(4,2,3,3)=(36*gmet(1,2)**3*dgmetds(1,3)+gmet(1,2)*(36*gmet(1,3)&
&       *gmet(2,2)*dgmetds(1,1)+gmet(1,1)*(288*gmet(2,3)*dgmetds(1,2)&
&       +180*gmet(2,2)*dgmetds(1,3)))+gmet(1,1)*gmet(2,2)*(144*gmet(2,3)&
&       *dgmetds(1,1)-108*gmet(1,3)*dgmetds(1,2)-36*gmet(1,1)*dgmetds(2,3))&
&       +gmet(1,2)**2*(36*gmet(2,3)*dgmetds(1,1)+36*gmet(1,3)*dgmetds(1,2)&
&       +108*gmet(1,1)*dgmetds(2,3)))/12.d0
       cm(5,2,3,3)=12*gmet(1,1)*gmet(1,2)*gmet(2,3)*dgmetds(1,1)+gmet(1,1)&
&       *gmet(2,2)*(-4.5d0*gmet(1,3)*dgmetds(1,1)-3*gmet(1,1)*dgmetds(1,3))&
&       +gmet(1,2)**2*(1.5d0*gmet(1,3)*dgmetds(1,1)+9*gmet(1,1)*dgmetds(1,3))
       cm(6,2,3,3)=1.5d0*gmet(1,2)**3*dgmetds(1,1)+7.5d0*gmet(1,1)*gmet(1,2)&
&       *gmet(2,2)*dgmetds(1,1)+9*gmet(1,1)*gmet(1,2)**2*dgmetds(1,2)&
&       -3*gmet(1,1)**2*gmet(2,2)*dgmetds(1,2)
       cm(7,2,3,3)=6*gmet(1,2)**2*gmet(2,2)*dgmetds(1,2)+12*gmet(1,1)&
&       *gmet(2,2)**2*dgmetds(1,2)+1.5d0*gmet(1,2)**3*dgmetds(2,2)+gmet(1,2)&
&       *gmet(2,2)*(3*gmet(2,2)*dgmetds(1,1)+7.5d0*gmet(1,1)*dgmetds(2,2))
       cm(8,2,3,3)=((48*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(6*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+2*(-12*gmet(1,3)**2*gmet(2,2)&
&       +36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)&
&       *(30*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+4*(6*gmet(1,2)&
&       *gmet(1,3)*gmet(2,2)+6*gmet(1,2)**2*gmet(2,3)+24*gmet(1,1)*gmet(2,2)&
&       *gmet(2,3))*dgmetds(1,3)+2*(6*gmet(1,2)**2*gmet(1,3)-18*gmet(1,1)&
&       *gmet(1,3)*gmet(2,2)+48*gmet(1,1)*gmet(1,2)*gmet(2,3))*dgmetds(2,3)&
&       +(6*gmet(1,2)**3+30*gmet(1,1)*gmet(1,2)*gmet(2,2))*dgmetds(3,3))&
&       /4.d0
       cm(9,2,3,3)=(gmet(2,2)*(24*gmet(1,3)*gmet(2,2)+12*gmet(1,2)*gmet(2,3))&
&       *dgmetds(1,1)+4*(6*gmet(1,2)*gmet(1,3)*gmet(2,2)+6*gmet(1,2)&
&       **2*gmet(2,3)+24*gmet(1,1)*gmet(2,2)*gmet(2,3))*dgmetds(1,2)&
&       +2*gmet(2,2)*(12*gmet(1,2)**2+24*gmet(1,1)*gmet(2,2))*dgmetds(1,3)&
&       +(6*gmet(1,2)**2*gmet(1,3)-18*gmet(1,1)*gmet(1,3)*gmet(2,2)+48*gmet(1,1)&
&       *gmet(1,2)*gmet(2,3))*dgmetds(2,2)+2*(6*gmet(1,2)**3+30*gmet(1,1)&
&       *gmet(1,2)*gmet(2,2))*dgmetds(2,3))/4.d0
       cm(10,2,3,3)=(2*(-36*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+12*(-12*gmet(1,3)**2*gmet(2,2)&
&       +36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)&
&       *(30*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,3)+6*(6*gmet(1,2)&
&       **2*gmet(1,3)-18*gmet(1,1)*gmet(1,3)*gmet(2,2)+48*gmet(1,1)*gmet(1,2)&
&       *gmet(2,3))*dgmetds(3,3))/24.d0
       cm(11,2,3,3)=gmet(2,2)*(6*gmet(1,2)*gmet(2,2)*dgmetds(1,2)+3*gmet(1,2)&
&       **2*dgmetds(2,2)+6*gmet(1,1)*gmet(2,2)*dgmetds(2,2))
       cm(12,2,3,3)=(2*(48*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(6*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+2*gmet(2,2)*(24*gmet(1,3)&
&       *gmet(2,2)+12*gmet(1,2)*gmet(2,3))*dgmetds(1,3)+(-12*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(30*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(2,2)&
&       +4*(6*gmet(1,2)*gmet(1,3)*gmet(2,2)+6*gmet(1,2)**2*gmet(2,3)&
&       +24*gmet(1,1)*gmet(2,2)*gmet(2,3))*dgmetds(2,3)+gmet(2,2)*(12*gmet(1,2)&
&       **2+24*gmet(1,1)*gmet(2,2))*dgmetds(3,3))/4.d0
       cm(13,2,3,3)=(72*gmet(1,2)*gmet(2,2)*(gmet(2,3)*dgmetds(1,2)+gmet(2,2)&
&       *dgmetds(1,3))+gmet(1,3)*gmet(2,2)*(144*gmet(2,2)*dgmetds(1,2)&
&       +36*gmet(1,2)*dgmetds(2,2))+144*gmet(1,1)*gmet(2,2)*(gmet(2,3)&
&       *dgmetds(2,2)+gmet(2,2)*dgmetds(2,3))+gmet(1,2)**2*(36*gmet(2,3)&
&       *dgmetds(2,2)+72*gmet(2,2)*dgmetds(2,3)))/12.d0
       cm(14,2,3,3)=(2*(-36*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+6*(48*gmet(1,3)*gmet(2,2)&
&       *gmet(2,3)+gmet(1,2)*(6*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))&
&       *dgmetds(1,3)+6*(-12*gmet(1,3)**2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(30*gmet(2,3)&
&       **2-6*gmet(2,2)*gmet(3,3)))*dgmetds(2,3)+6*(6*gmet(1,2)*gmet(1,3)&
&       *gmet(2,2)+6*gmet(1,2)**2*gmet(2,3)+24*gmet(1,1)*gmet(2,2)*gmet(2,3))&
&       *dgmetds(3,3))/12.d0
       cm(15,2,3,3)=(8*(-36*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,3)+12*(-12*gmet(1,3)**2*gmet(2,2)&
&       +36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)&
&       *(30*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(3,3))/48.d0
       cm(16,2,3,3)=3*gmet(1,2)*gmet(2,2)**2*dgmetds(2,2)
       cm(17,2,3,3)=(288*gmet(1,3)*gmet(2,2)*(gmet(2,3)*dgmetds(2,2)&
&       +gmet(2,2)*dgmetds(2,3))+gmet(1,2)*(36*gmet(2,3)**2*dgmetds(2,2)&
&       +144*gmet(2,2)*gmet(2,3)*dgmetds(2,3)+gmet(2,2)*(-108*gmet(3,3)&
&       *dgmetds(2,2)+72*gmet(2,2)*dgmetds(3,3))))/24.d0
       cm(18,2,3,3)=gmet(2,2)*(6*gmet(1,3)*gmet(2,2)*dgmetds(2,2)+gmet(1,2)&
&       *(3*gmet(2,3)*dgmetds(2,2)+6*gmet(2,2)*dgmetds(2,3)))
       cm(19,2,3,3)=(2*(-36*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(2,2)+12*(48*gmet(1,3)*gmet(2,2)&
&       *gmet(2,3)+gmet(1,2)*(6*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))&
&       *dgmetds(2,3)+6*gmet(2,2)*(24*gmet(1,3)*gmet(2,2)+12*gmet(1,2)&
&       *gmet(2,3))*dgmetds(3,3))/24.0d0
       cm(20,2,3,3)=(8*(-36*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(2,3)+12*(48*gmet(1,3)*gmet(2,2)&
&       *gmet(2,3)+gmet(1,2)*(6*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))&
&       *dgmetds(3,3))/48.d0
       cm(21,2,3,3)=((-36*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(3,3))/12.d0
       cm(1,3,3,3)=gmet(1,1)*(4.5d0*gmet(1,3)**2-1.5d0*gmet(1,1)*gmet(3,3))&
&       *dgmetds(1,1)
       cm(2,3,3,3)=(6*(-12*gmet(1,3)**2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(30*gmet(2,3)&
&       **2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+12*(48*gmet(1,1)*gmet(1,3)&
&       *gmet(2,3)+gmet(1,2)*(6*gmet(1,3)**2-18*gmet(1,1)*gmet(3,3)))&
&       *dgmetds(1,2)+2*gmet(1,1)*(54*gmet(1,3)**2-18*gmet(1,1)*gmet(3,3))&
&       *dgmetds(2,2))/24.d0
       cm(3,3,3,3)=3*gmet(1,3)**3*dgmetds(1,3)+15*gmet(1,1)*gmet(1,3)&
&       *gmet(3,3)*dgmetds(1,3)+gmet(1,1)*gmet(3,3)*(6*gmet(3,3)*dgmetds(1,1)&
&       -1.5d0*gmet(1,1)*dgmetds(3,3))+gmet(1,3)**2*(3*gmet(3,3)*dgmetds(1,1)&
&       +4.5d0*gmet(1,1)*dgmetds(3,3))
       cm(4,3,3,3)=(36*gmet(1,3)**3*dgmetds(1,2)+gmet(1,3)*(36*gmet(1,2)&
&       *gmet(3,3)*dgmetds(1,1)+gmet(1,1)*(180*gmet(3,3)*dgmetds(1,2)&
&       +288*gmet(2,3)*dgmetds(1,3)))+gmet(1,1)*gmet(3,3)*(144*gmet(2,3)&
&       *dgmetds(1,1)-108*gmet(1,2)*dgmetds(1,3)-36*gmet(1,1)*dgmetds(2,3))&
&       +gmet(1,3)**2*(36*gmet(2,3)*dgmetds(1,1)+36*gmet(1,2)*dgmetds(1,3)&
&       +108*gmet(1,1)*dgmetds(2,3)))/12.d0
       cm(5,3,3,3)=1.5d0*gmet(1,3)**3*dgmetds(1,1)+7.5d0*gmet(1,1)*gmet(1,3)&
&       *gmet(3,3)*dgmetds(1,1)+9*gmet(1,1)*gmet(1,3)**2*dgmetds(1,3)&
&       -3*gmet(1,1)**2*gmet(3,3)*dgmetds(1,3)
       cm(6,3,3,3)=gmet(1,2)*(1.5d0*gmet(1,3)**2-4.5d0*gmet(1,1)*gmet(3,3))&
&       *dgmetds(1,1)+gmet(1,1)*(12*gmet(1,3)*gmet(2,3)*dgmetds(1,1)&
&       +9*gmet(1,3)**2*dgmetds(1,2)-3*gmet(1,1)*gmet(3,3)*dgmetds(1,2))
       cm(7,3,3,3)=(2*(-36*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+12*(-12*gmet(1,3)**2*gmet(2,2)&
&       +36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)&
&       *(30*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+6*(48*gmet(1,1)&
&       *gmet(1,3)*gmet(2,3)+gmet(1,2)*(6*gmet(1,3)**2-18*gmet(1,1)*gmet(3,3)))&
&       *dgmetds(2,2))/24.d0
       cm(8,3,3,3)=(gmet(3,3)*(12*gmet(1,3)*gmet(2,3)+24*gmet(1,2)*gmet(3,3))&
&       *dgmetds(1,1)+2*gmet(3,3)*(12*gmet(1,3)**2+24*gmet(1,1)*gmet(3,3))&
&       *dgmetds(1,2)+4*(6*gmet(1,3)**2*gmet(2,3)+6*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)+24*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(1,3)+2*(6*gmet(1,3)&
&       **3+30*gmet(1,1)*gmet(1,3)*gmet(3,3))*dgmetds(2,3)+(48*gmet(1,1)&
&       *gmet(1,3)*gmet(2,3)+gmet(1,2)*(6*gmet(1,3)**2-18*gmet(1,1)*gmet(3,3)))&
&       *dgmetds(3,3))/4.d0
       cm(9,3,3,3)=((48*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(6*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+4*(6*gmet(1,3)**2*gmet(2,3)&
&       +6*gmet(1,2)*gmet(1,3)*gmet(3,3)+24*gmet(1,1)*gmet(2,3)*gmet(3,3))&
&       *dgmetds(1,2)+2*(-12*gmet(1,3)**2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(30*gmet(2,3)&
&       **2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,3)+(6*gmet(1,3)**3+30*gmet(1,1)&
&       *gmet(1,3)*gmet(3,3))*dgmetds(2,2)+2*(48*gmet(1,1)*gmet(1,3)&
&       *gmet(2,3)+gmet(1,2)*(6*gmet(1,3)**2-18*gmet(1,1)*gmet(3,3)))&
&       *dgmetds(2,3))/4.d0
       cm(10,3,3,3)=6*gmet(1,3)**2*gmet(3,3)*dgmetds(1,3)+12*gmet(1,1)&
&       *gmet(3,3)**2*dgmetds(1,3)+1.5d0*gmet(1,3)**3*dgmetds(3,3)+gmet(1,3)&
&       *gmet(3,3)*(3*gmet(3,3)*dgmetds(1,1)+7.5d0*gmet(1,1)*dgmetds(3,3))
       cm(11,3,3,3)=(8*(-36*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+12*(-12*gmet(1,3)**2*gmet(2,2)&
&       +36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)&
&       *(30*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(2,2))/48.d0
       cm(12,3,3,3)=(2*gmet(3,3)*(12*gmet(1,3)*gmet(2,3)+24*gmet(1,2)&
&       *gmet(3,3))*dgmetds(1,2)+2*(48*gmet(1,2)*gmet(2,3)*gmet(3,3)&
&       +gmet(1,3)*(6*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,3)&
&       +gmet(3,3)*(12*gmet(1,3)**2+24*gmet(1,1)*gmet(3,3))*dgmetds(2,2)&
&       +4*(6*gmet(1,3)**2*gmet(2,3)+6*gmet(1,2)*gmet(1,3)*gmet(3,3)&
&       +24*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(2,3)+(-12*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(30*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(3,3))&
&       /4.d0
       cm(13,3,3,3)=(6*(48*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(6*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+2*(-36*gmet(1,3)*gmet(2,2)&
&       *gmet(2,3)+gmet(1,2)*(90*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))&
&       *dgmetds(1,3)+6*(6*gmet(1,3)**2*gmet(2,3)+6*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)+24*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(2,2)+6*(-12*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(30*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(2,3))&
&       /12.d0
       cm(14,3,3,3)=(gmet(1,3)*gmet(3,3)*(72*gmet(3,3)*dgmetds(1,2)+72*gmet(2,3)&
&       *dgmetds(1,3)+36*gmet(1,2)*dgmetds(3,3))+gmet(1,3)**2*(72*gmet(3,3)&
&       *dgmetds(2,3)+36*gmet(2,3)*dgmetds(3,3))+144*gmet(3,3)*(gmet(1,2)&
&       *gmet(3,3)*dgmetds(1,3)+gmet(1,1)*(gmet(3,3)*dgmetds(2,3)+gmet(2,3)&
&       *dgmetds(3,3))))/12.d0
       cm(15,3,3,3)=gmet(3,3)*(6*gmet(1,3)*gmet(3,3)*dgmetds(1,3)+3*gmet(1,3)&
&       **2*dgmetds(3,3)+6*gmet(1,1)*gmet(3,3)*dgmetds(3,3))
       cm(16,3,3,3)=((-36*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(2,2))/12.d0
       cm(17,3,3,3)=(6*gmet(3,3)*(12*gmet(1,3)*gmet(2,3)+24*gmet(1,2)&
&       *gmet(3,3))*dgmetds(2,2)+12*(48*gmet(1,2)*gmet(2,3)*gmet(3,3)&
&       +gmet(1,3)*(6*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))*dgmetds(2,3)&
&       +2*(-36*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(3,3))/24.d0
       cm(18,3,3,3)=(12*(48*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(6*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(2,2)+8*(-36*gmet(1,3)*gmet(2,2)&
&       *gmet(2,3)+gmet(1,2)*(90*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))&
&       *dgmetds(2,3))/48.d0
       cm(19,3,3,3)=(288*gmet(1,2)*gmet(3,3)*(gmet(3,3)*dgmetds(2,3)&
&       +gmet(2,3)*dgmetds(3,3))+gmet(1,3)*(72*gmet(3,3)**2*dgmetds(2,2)&
&       +36*gmet(2,3)**2*dgmetds(3,3)+gmet(3,3)*(144*gmet(2,3)*dgmetds(2,3)&
&       -108*gmet(2,2)*dgmetds(3,3))))/24.d0
       cm(20,3,3,3)=gmet(3,3)*(6*gmet(1,2)*gmet(3,3)*dgmetds(3,3)+gmet(1,3)&
&       *(6*gmet(3,3)*dgmetds(2,3)+3*gmet(2,3)*dgmetds(3,3)))
       cm(21,3,3,3)=3*gmet(1,3)*gmet(3,3)**2*dgmetds(3,3)
       cm(1,4,3,3)=(gmet(1,1)*(54*gmet(1,2)*gmet(1,3)-18*gmet(1,1)*gmet(2,3))&
&       *dgmetds(1,1))/6.d0
       cm(2,4,3,3)=(6*(6*gmet(1,2)*gmet(1,3)*gmet(2,2)+6*gmet(1,2)**2*gmet(2,3)&
&       +24*gmet(1,1)*gmet(2,2)*gmet(2,3))*dgmetds(1,1)+12*(6*gmet(1,2)&
&       **2*gmet(1,3)+24*gmet(1,1)*gmet(1,3)*gmet(2,2)+6*gmet(1,1)*gmet(1,2)&
&       *gmet(2,3))*dgmetds(1,2)+2*gmet(1,1)*(54*gmet(1,2)*gmet(1,3)&
&       -18*gmet(1,1)*gmet(2,3))*dgmetds(2,2))/12.0d0
       cm(3,4,3,3)=(6*(6*gmet(1,3)**2*gmet(2,3)+6*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)+24*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(1,1)+12*(6*gmet(1,1)&
&       *gmet(1,3)*gmet(2,3)+gmet(1,2)*(6*gmet(1,3)**2+24*gmet(1,1)*gmet(3,3)))&
&       *dgmetds(1,3)+2*gmet(1,1)*(54*gmet(1,2)*gmet(1,3)-18*gmet(1,1)&
&       *gmet(2,3))*dgmetds(3,3))/12.0d0
       cm(4,4,3,3)=(6*(9*gmet(1,3)**2*gmet(2,2)-6*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)+9*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(9*gmet(2,3)**2+15*gmet(2,2)&
&       *gmet(3,3)))*dgmetds(1,1)+6*(6*gmet(1,1)*gmet(1,3)*gmet(2,3)&
&       +gmet(1,2)*(6*gmet(1,3)**2+24*gmet(1,1)*gmet(3,3)))*dgmetds(1,2)&
&       +6*(6*gmet(1,2)**2*gmet(1,3)+24*gmet(1,1)*gmet(1,3)*gmet(2,2)&
&       +6*gmet(1,1)*gmet(1,2)*gmet(2,3))*dgmetds(1,3)+2*gmet(1,1)*(54*gmet(1,2)&
&       *gmet(1,3)-18*gmet(1,1)*gmet(2,3))*dgmetds(2,3))/6.0d0
       cm(5,4,3,3)=gmet(1,1)*gmet(2,3)*(3*gmet(1,3)*dgmetds(1,1)-6*gmet(1,1)&
&       *dgmetds(1,3))+gmet(1,2)*(3*gmet(1,3)**2*dgmetds(1,1)+12*gmet(1,1)&
&       *gmet(3,3)*dgmetds(1,1)+18*gmet(1,1)*gmet(1,3)*dgmetds(1,3))
       cm(6,4,3,3)=3*gmet(1,2)**2*gmet(1,3)*dgmetds(1,1)+gmet(1,1)*gmet(1,2)&
&       *(3*gmet(2,3)*dgmetds(1,1)+18*gmet(1,3)*dgmetds(1,2))+gmet(1,1)&
&       *(12*gmet(1,3)*gmet(2,2)*dgmetds(1,1)-6*gmet(1,1)*gmet(2,3)*dgmetds(1,2))
       cm(7,4,3,3)=(2*gmet(2,2)*(-18*gmet(1,3)*gmet(2,2)+54*gmet(1,2)&
&       *gmet(2,3))*dgmetds(1,1)+12*(6*gmet(1,2)*gmet(1,3)*gmet(2,2)&
&       +6*gmet(1,2)**2*gmet(2,3)+24*gmet(1,1)*gmet(2,2)*gmet(2,3))*dgmetds(1,2)&
&       +6*(6*gmet(1,2)**2*gmet(1,3)+24*gmet(1,1)*gmet(1,3)*gmet(2,2)&
&       +6*gmet(1,1)*gmet(1,2)*gmet(2,3))*dgmetds(2,2))/12.d0
       cm(8,4,3,3)=((6*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(6*gmet(2,3)&
&       **2+24*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+2*(6*gmet(1,3)**2*gmet(2,3)&
&       +6*gmet(1,2)*gmet(1,3)*gmet(3,3)+24*gmet(1,1)*gmet(2,3)*gmet(3,3))&
&       *dgmetds(1,2)+4*(9*gmet(1,3)**2*gmet(2,2)-6*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)+9*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(9*gmet(2,3)**2+15*gmet(2,2)&
&       *gmet(3,3)))*dgmetds(1,3)+2*(6*gmet(1,1)*gmet(1,3)*gmet(2,3)&
&       +gmet(1,2)*(6*gmet(1,3)**2+24*gmet(1,1)*gmet(3,3)))*dgmetds(2,3)&
&       +(6*gmet(1,2)**2*gmet(1,3)+24*gmet(1,1)*gmet(1,3)*gmet(2,2)+6*gmet(1,1)&
&       *gmet(1,2)*gmet(2,3))*dgmetds(3,3))/2.d0
       cm(9,4,3,3)=((6*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(6*gmet(2,3)&
&       **2+24*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+4*(9*gmet(1,3)**2*gmet(2,2)&
&       -6*gmet(1,2)*gmet(1,3)*gmet(2,3)+9*gmet(1,2)**2*gmet(3,3)+gmet(1,1)&
&       *(9*gmet(2,3)**2+15*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+2*(6*gmet(1,2)&
&       *gmet(1,3)*gmet(2,2)+6*gmet(1,2)**2*gmet(2,3)+24*gmet(1,1)*gmet(2,2)&
&       *gmet(2,3))*dgmetds(1,3)+(6*gmet(1,1)*gmet(1,3)*gmet(2,3)+gmet(1,2)&
&       *(6*gmet(1,3)**2+24*gmet(1,1)*gmet(3,3)))*dgmetds(2,2)+2*(6*gmet(1,2)&
&       **2*gmet(1,3)+24*gmet(1,1)*gmet(1,3)*gmet(2,2)+6*gmet(1,1)*gmet(1,2)&
&       *gmet(2,3))*dgmetds(2,3))/2.d0
       cm(10,4,3,3)=(2*gmet(3,3)*(54*gmet(1,3)*gmet(2,3)-18*gmet(1,2)&
&       *gmet(3,3))*dgmetds(1,1)+12*(6*gmet(1,3)**2*gmet(2,3)+6*gmet(1,2)&
&       *gmet(1,3)*gmet(3,3)+24*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(1,3)&
&       +6*(6*gmet(1,1)*gmet(1,3)*gmet(2,3)+gmet(1,2)*(6*gmet(1,3)**2+24*gmet(1,1)&
&       *gmet(3,3)))*dgmetds(3,3))/12.d0
       cm(11,4,3,3)=gmet(1,3)*gmet(2,2)*(-6*gmet(2,2)*dgmetds(1,2)+3*gmet(1,2)&
&       *dgmetds(2,2))+gmet(2,3)*(18*gmet(1,2)*gmet(2,2)*dgmetds(1,2)&
&       +3*gmet(1,2)**2*dgmetds(2,2)+12*gmet(1,1)*gmet(2,2)*dgmetds(2,2))
       cm(12,4,3,3)=(2*(6*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(6*gmet(2,3)&
&       **2+24*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+2*(6*gmet(1,3)*gmet(2,2)&
&       *gmet(2,3)+gmet(1,2)*(6*gmet(2,3)**2+24*gmet(2,2)*gmet(3,3)))&
&       *dgmetds(1,3)+(6*gmet(1,3)**2*gmet(2,3)+6*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)+24*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(2,2)+4*(9*gmet(1,3)&
&       **2*gmet(2,2)-6*gmet(1,2)*gmet(1,3)*gmet(2,3)+9*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(9*gmet(2,3)**2+15*gmet(2,2)*gmet(3,3)))*dgmetds(2,3)&
&       +(6*gmet(1,2)*gmet(1,3)*gmet(2,2)+6*gmet(1,2)**2*gmet(2,3)+24*gmet(1,1)&
&       *gmet(2,2)*gmet(2,3))*dgmetds(3,3))/2.d0
       cm(13,4,3,3)=(6*(6*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(6*gmet(2,3)&
&       **2+24*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+2*gmet(2,2)*(-18*gmet(1,3)&
&       *gmet(2,2)+54*gmet(1,2)*gmet(2,3))*dgmetds(1,3)+6*(9*gmet(1,3)&
&       **2*gmet(2,2)-6*gmet(1,2)*gmet(1,3)*gmet(2,3)+9*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(9*gmet(2,3)**2+15*gmet(2,2)*gmet(3,3)))*dgmetds(2,2)&
&       +6*(6*gmet(1,2)*gmet(1,3)*gmet(2,2)+6*gmet(1,2)**2*gmet(2,3)&
&       +24*gmet(1,1)*gmet(2,2)*gmet(2,3))*dgmetds(2,3))/6.d0
       cm(14,4,3,3)=(2*gmet(3,3)*(54*gmet(1,3)*gmet(2,3)-18*gmet(1,2)&
&       *gmet(3,3))*dgmetds(1,2)+6*(6*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)&
&       *(6*gmet(2,3)**2+24*gmet(2,2)*gmet(3,3)))*dgmetds(1,3)+6*(6*gmet(1,3)&
&       **2*gmet(2,3)+6*gmet(1,2)*gmet(1,3)*gmet(3,3)+24*gmet(1,1)*gmet(2,3)&
&       *gmet(3,3))*dgmetds(2,3)+6*(9*gmet(1,3)**2*gmet(2,2)-6*gmet(1,2)&
&       *gmet(1,3)*gmet(2,3)+9*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(9*gmet(2,3)&
&       **2+15*gmet(2,2)*gmet(3,3)))*dgmetds(3,3))/6.d0
       cm(15,4,3,3)=3*gmet(1,3)**2*gmet(2,3)*dgmetds(3,3)+gmet(1,3)*gmet(3,3)&
&       *(18*gmet(2,3)*dgmetds(1,3)+3*gmet(1,2)*dgmetds(3,3))+gmet(3,3)&
&       *(-6*gmet(1,2)*gmet(3,3)*dgmetds(1,3)+12*gmet(1,1)*gmet(2,3)&
&       *dgmetds(3,3))
       cm(16,4,3,3)=(gmet(2,2)*(-18*gmet(1,3)*gmet(2,2)+54*gmet(1,2)&
&       *gmet(2,3))*dgmetds(2,2))/6.d0
       cm(17,4,3,3)=(6*(6*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(6*gmet(2,3)&
&       **2+24*gmet(2,2)*gmet(3,3)))*dgmetds(2,2)+12*(6*gmet(1,3)*gmet(2,2)&
&       *gmet(2,3)+gmet(1,2)*(6*gmet(2,3)**2+24*gmet(2,2)*gmet(3,3)))&
&       *dgmetds(2,3)+2*gmet(2,2)*(-18*gmet(1,3)*gmet(2,2)+54*gmet(1,2)&
&       *gmet(2,3))*dgmetds(3,3))/12.0d0
       cm(18,4,3,3)=gmet(1,3)*gmet(2,2)*(3*gmet(2,3)*dgmetds(2,2)-6*gmet(2,2)&
&       *dgmetds(2,3))+gmet(1,2)*(3*gmet(2,3)**2*dgmetds(2,2)+12*gmet(2,2)&
&       *gmet(3,3)*dgmetds(2,2)+18*gmet(2,2)*gmet(2,3)*dgmetds(2,3))
       cm(19,4,3,3)=(2*gmet(3,3)*(54*gmet(1,3)*gmet(2,3)-18*gmet(1,2)&
&       *gmet(3,3))*dgmetds(2,2)+12*(6*gmet(1,2)*gmet(2,3)*gmet(3,3)&
&       +gmet(1,3)*(6*gmet(2,3)**2+24*gmet(2,2)*gmet(3,3)))*dgmetds(2,3)&
&       +6*(6*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(6*gmet(2,3)**2+24*gmet(2,2)&
&       *gmet(3,3)))*dgmetds(3,3))/12.d0
       cm(20,4,3,3)=gmet(1,2)*gmet(3,3)*(-6*gmet(3,3)*dgmetds(2,3)+3*gmet(2,3)&
&       *dgmetds(3,3))+gmet(1,3)*(18*gmet(2,3)*gmet(3,3)*dgmetds(2,3)&
&       +3*gmet(2,3)**2*dgmetds(3,3)+12*gmet(2,2)*gmet(3,3)*dgmetds(3,3))
       cm(21,4,3,3)=(gmet(3,3)*(54*gmet(1,3)*gmet(2,3)-18*gmet(1,2)*gmet(3,3))&
&       *dgmetds(3,3))/6.d0
       cm(1,5,3,3)=3*gmet(1,1)**2*gmet(1,3)*dgmetds(1,1)
       cm(2,5,3,3)=1.5d0*gmet(1,2)**2*gmet(1,3)*dgmetds(1,1)+gmet(1,1)&
&       *gmet(1,2)*(12*gmet(2,3)*dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,2))&
&       +gmet(1,1)*(12*gmet(1,1)*gmet(2,3)*dgmetds(1,2)+gmet(1,3)*(-4.5d0*gmet(2,2)&
&       *dgmetds(1,1)+3*gmet(1,1)*dgmetds(2,2)))
       cm(3,5,3,3)=1.5d0*gmet(1,3)**3*dgmetds(1,1)+6*gmet(1,1)*gmet(1,3)&
&       **2*dgmetds(1,3)+12*gmet(1,1)**2*gmet(3,3)*dgmetds(1,3)+gmet(1,1)&
&       *gmet(1,3)*(7.5d0*gmet(3,3)*dgmetds(1,1)+3*gmet(1,1)*dgmetds(3,3))
       cm(4,5,3,3)=gmet(1,2)*(3*gmet(1,3)**2*dgmetds(1,1)+12*gmet(1,1)&
&       *gmet(3,3)*dgmetds(1,1)+6*gmet(1,1)*gmet(1,3)*dgmetds(1,3))+gmet(1,1)&
&       *(6*gmet(1,3)**2*dgmetds(1,2)+12*gmet(1,1)*(gmet(3,3)*dgmetds(1,2)&
&       +gmet(2,3)*dgmetds(1,3))+gmet(1,3)*(3*gmet(2,3)*dgmetds(1,1)&
&       +6*gmet(1,1)*dgmetds(2,3)))
       cm(5,5,3,3)=gmet(1,1)*(3*gmet(1,3)**2*dgmetds(1,1)+6*gmet(1,1)&
&       *gmet(3,3)*dgmetds(1,1)+6*gmet(1,1)*gmet(1,3)*dgmetds(1,3))
       cm(6,5,3,3)=gmet(1,1)*(3*gmet(1,2)*gmet(1,3)*dgmetds(1,1)+6*gmet(1,1)&
&       *(gmet(2,3)*dgmetds(1,1)+gmet(1,3)*dgmetds(1,2)))
       cm(7,5,3,3)=gmet(1,2)**2*(7.5d0*gmet(2,3)*dgmetds(1,1)+3*gmet(1,3)&
&       *dgmetds(1,2))+gmet(1,1)*(gmet(2,2)*(-1.5d0*gmet(2,3)*dgmetds(1,1)&
&       -9*gmet(1,3)*dgmetds(1,2))+6*gmet(1,1)*gmet(2,3)*dgmetds(2,2))&
&       +gmet(1,2)*(24*gmet(1,1)*gmet(2,3)*dgmetds(1,2)+gmet(1,3)*(-3*gmet(2,2)&
&       *dgmetds(1,1)+3*gmet(1,1)*dgmetds(2,2)))
       cm(8,5,3,3)=3*gmet(1,3)**3*dgmetds(1,2)+gmet(1,3)**2*(1.5d0*gmet(2,3)&
&       *dgmetds(1,1)+6*(gmet(1,2)*dgmetds(1,3)+gmet(1,1)*dgmetds(2,3)))&
&       +gmet(1,3)*(gmet(1,1)*(15*gmet(3,3)*dgmetds(1,2)+6*gmet(2,3)&
&       *dgmetds(1,3))+gmet(1,2)*(12*gmet(3,3)*dgmetds(1,1)+3*gmet(1,1)&
&       *dgmetds(3,3)))+gmet(1,1)*(gmet(3,3)*(24*gmet(1,2)*dgmetds(1,3)&
&       +12*gmet(1,1)*dgmetds(2,3))+gmet(2,3)*(-4.5d0*gmet(3,3)*dgmetds(1,1)&
&       +6*gmet(1,1)*dgmetds(3,3)))
       cm(9,5,3,3)=7.5d0*gmet(1,2)**2*gmet(3,3)*dgmetds(1,1)+gmet(1,1)&
&       *(-3*gmet(2,3)**2*dgmetds(1,1)-1.5d0*gmet(2,2)*gmet(3,3)*dgmetds(1,1)&
&       +24*gmet(1,2)*gmet(3,3)*dgmetds(1,2)+24*gmet(1,2)*gmet(2,3)*dgmetds(1,3))&
&       +gmet(1,3)**2*(-3*gmet(2,2)*dgmetds(1,1)+6*gmet(1,2)*dgmetds(1,2)&
&       +3*gmet(1,1)*dgmetds(2,2))+gmet(1,1)**2*(6*gmet(3,3)*dgmetds(2,2)&
&       +12*gmet(2,3)*dgmetds(2,3))+gmet(1,3)*(3*gmet(1,2)**2*dgmetds(1,3)&
&       +gmet(1,1)*(6*gmet(2,3)*dgmetds(1,2)-9*gmet(2,2)*dgmetds(1,3))&
&       +gmet(1,2)*(9*gmet(2,3)*dgmetds(1,1)+6*gmet(1,1)*dgmetds(2,3)))
       cm(10,5,3,3)=3*gmet(1,3)**3*dgmetds(1,3)+15*gmet(1,1)*gmet(1,3)&
&       *gmet(3,3)*dgmetds(1,3)+gmet(1,3)**2*(4.5d0*gmet(3,3)*dgmetds(1,1)&
&       +3*gmet(1,1)*dgmetds(3,3))+gmet(1,1)*gmet(3,3)*(-1.5d0*gmet(3,3)&
&       *dgmetds(1,1)+6*gmet(1,1)*dgmetds(3,3))
       cm(11,5,3,3)=gmet(1,1)*gmet(2,2)*(-3*gmet(2,3)*dgmetds(1,2)-4.5d0*gmet(1,3)&
&       *dgmetds(2,2))+gmet(1,2)**2*(15*gmet(2,3)*dgmetds(1,2)+1.5d0*gmet(1,3)&
&       *dgmetds(2,2))+gmet(1,2)*(-6*gmet(1,3)*gmet(2,2)*dgmetds(1,2)&
&       +12*gmet(1,1)*gmet(2,3)*dgmetds(2,2))
       cm(12,5,3,3)=15*gmet(1,2)**2*gmet(3,3)*dgmetds(1,3)+1.5d0*gmet(1,3)&
&       **3*dgmetds(2,2)+gmet(1,3)**2*(3*gmet(2,3)*dgmetds(1,2)-6*gmet(2,2)&
&       *dgmetds(1,3)+6*gmet(1,2)*dgmetds(2,3))+gmet(1,1)*(-9*gmet(2,3)&
&       *gmet(3,3)*dgmetds(1,2)-6*gmet(2,3)**2*dgmetds(1,3)-3*gmet(2,2)&
&       *gmet(3,3)*dgmetds(1,3)+24*gmet(1,2)*gmet(3,3)*dgmetds(2,3)+12*gmet(1,2)&
&       *gmet(2,3)*dgmetds(3,3))+gmet(1,3)*(gmet(1,2)*(24*gmet(3,3)*dgmetds(1,2)&
&       +18*gmet(2,3)*dgmetds(1,3))+1.5d0*gmet(1,2)**2*dgmetds(3,3)+gmet(1,1)&
&       *(7.5d0*gmet(3,3)*dgmetds(2,2)+6*gmet(2,3)*dgmetds(2,3)-4.5d0*gmet(2,2)&
&       *dgmetds(3,3)))
       cm(13,5,3,3)=15*gmet(1,2)**2*(gmet(3,3)*dgmetds(1,2)+gmet(2,3)&
&       *dgmetds(1,3))+gmet(1,3)**2*(-6*gmet(2,2)*dgmetds(1,2)+3*gmet(1,2)&
&       *dgmetds(2,2))+gmet(1,1)*(-6*gmet(2,3)**2*dgmetds(1,2)-3*gmet(2,2)&
&       *gmet(3,3)*dgmetds(1,2)-3*gmet(2,2)*gmet(2,3)*dgmetds(1,3)+12*gmet(1,2)&
&       *gmet(3,3)*dgmetds(2,2)+24*gmet(1,2)*gmet(2,3)*dgmetds(2,3))&
&       +gmet(1,3)*(gmet(1,2)*(18*gmet(2,3)*dgmetds(1,2)-6*gmet(2,2)&
&       *dgmetds(1,3))+3*gmet(1,2)**2*dgmetds(2,3)+gmet(1,1)*(3*gmet(2,3)&
&       *dgmetds(2,2)-9*gmet(2,2)*dgmetds(2,3)))
       cm(14,5,3,3)=3*gmet(1,3)**3*dgmetds(2,3)+gmet(1,1)*gmet(3,3)*(-3*gmet(3,3)&
&       *dgmetds(1,2)-9*gmet(2,3)*dgmetds(1,3)+12*gmet(1,2)*dgmetds(3,3))&
&       +gmet(1,3)**2*(9*gmet(3,3)*dgmetds(1,2)+3*(gmet(2,3)*dgmetds(1,3)&
&       +gmet(1,2)*dgmetds(3,3)))+gmet(1,3)*(24*gmet(1,2)*gmet(3,3)*dgmetds(1,3)&
&       +gmet(1,1)*(15*gmet(3,3)*dgmetds(2,3)+3*gmet(2,3)*dgmetds(3,3)))
       cm(15,5,3,3)=9*gmet(1,3)**2*gmet(3,3)*dgmetds(1,3)-3*gmet(1,1)&
&       *gmet(3,3)**2*dgmetds(1,3)+1.5d0*gmet(1,3)**3*dgmetds(3,3)+7.5d0*gmet(1,1)&
&       *gmet(1,3)*gmet(3,3)*dgmetds(3,3)
       cm(16,5,3,3)=((-36*gmet(1,2)*gmet(1,3)*gmet(2,2)+90*gmet(1,2)&
&       **2*gmet(2,3)-18*gmet(1,1)*gmet(2,2)*gmet(2,3))*dgmetds(2,2))&
&       /12.d0
       cm(17,5,3,3)=gmet(1,3)**2*(1.5d0*gmet(2,3)*dgmetds(2,2)-6*gmet(2,2)&
&       *dgmetds(2,3))+gmet(1,2)*gmet(1,3)*(12*gmet(3,3)*dgmetds(2,2)&
&       +18*gmet(2,3)*dgmetds(2,3)-3*gmet(2,2)*dgmetds(3,3))+gmet(1,2)&
&       **2*(15*gmet(3,3)*dgmetds(2,3)+7.5d0*gmet(2,3)*dgmetds(3,3))&
&       +gmet(1,1)*(-4.5d0*gmet(2,3)*gmet(3,3)*dgmetds(2,2)-6*gmet(2,3)&
&       **2*dgmetds(2,3)-3*gmet(2,2)*gmet(3,3)*dgmetds(2,3)-1.5d0*gmet(2,2)&
&       *gmet(2,3)*dgmetds(3,3))
       cm(18,5,3,3)=(12*(-12*gmet(1,3)**2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)+30*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(-12*gmet(2,3)&
&       **2-6*gmet(2,2)*gmet(3,3)))*dgmetds(2,2)+8*(-36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,2)+90*gmet(1,2)**2*gmet(2,3)-18*gmet(1,1)*gmet(2,2)*gmet(2,3))&
&       *dgmetds(2,3))/48.d0
       cm(19,5,3,3)=7.5d0*gmet(1,2)**2*gmet(3,3)*dgmetds(3,3)+gmet(1,3)&
&       **2*(4.5d0*gmet(3,3)*dgmetds(2,2)+3*gmet(2,3)*dgmetds(2,3)-3*gmet(2,2)&
&       *dgmetds(3,3))+gmet(1,2)*gmet(1,3)*(24*gmet(3,3)*dgmetds(2,3)&
&       +9*gmet(2,3)*dgmetds(3,3))+gmet(1,1)*(-1.5d0*gmet(3,3)**2*dgmetds(2,2)&
&       -9*gmet(2,3)*gmet(3,3)*dgmetds(2,3)-3*gmet(2,3)**2*dgmetds(3,3)&
&       -1.5d0*gmet(2,2)*gmet(3,3)*dgmetds(3,3))
       cm(20,5,3,3)=12*gmet(1,2)*gmet(1,3)*gmet(3,3)*dgmetds(3,3)+gmet(1,1)&
&       *gmet(3,3)*(-3*gmet(3,3)*dgmetds(2,3)-4.5d0*gmet(2,3)*dgmetds(3,3))&
&       +gmet(1,3)**2*(9*gmet(3,3)*dgmetds(2,3)+1.5d0*gmet(2,3)*dgmetds(3,3))
       cm(21,5,3,3)=gmet(3,3)*(4.5d0*gmet(1,3)**2-1.5d0*gmet(1,1)*gmet(3,3))&
&       *dgmetds(3,3)
       cm(1,6,3,3)=3*gmet(1,1)**2*gmet(1,2)*dgmetds(1,1)
       cm(2,6,3,3)=1.5d0*gmet(1,2)**3*dgmetds(1,1)+6*gmet(1,1)*gmet(1,2)&
&       **2*dgmetds(1,2)+12*gmet(1,1)**2*gmet(2,2)*dgmetds(1,2)+gmet(1,1)&
&       *gmet(1,2)*(7.5d0*gmet(2,2)*dgmetds(1,1)+3*gmet(1,1)*dgmetds(2,2))
       cm(3,6,3,3)=12*gmet(1,1)*gmet(2,3)*(gmet(1,3)*dgmetds(1,1)+gmet(1,1)&
&       *dgmetds(1,3))+gmet(1,2)*(1.5d0*gmet(1,3)**2*dgmetds(1,1)+6*gmet(1,1)&
&       *gmet(1,3)*dgmetds(1,3)+gmet(1,1)*(-4.5d0*gmet(3,3)*dgmetds(1,1)&
&       +3*gmet(1,1)*dgmetds(3,3)))
       cm(4,6,3,3)=gmet(1,2)**2*(3*gmet(1,3)*dgmetds(1,1)+6*gmet(1,1)&
&       *dgmetds(1,3))+12*gmet(1,1)*(gmet(1,3)*gmet(2,2)*dgmetds(1,1)&
&       +gmet(1,1)*(gmet(2,3)*dgmetds(1,2)+gmet(2,2)*dgmetds(1,3)))+gmet(1,1)&
&       *gmet(1,2)*(3*gmet(2,3)*dgmetds(1,1)+6*(gmet(1,3)*dgmetds(1,2)&
&       +gmet(1,1)*dgmetds(2,3)))
       cm(5,6,3,3)=gmet(1,1)*(6*gmet(1,1)*gmet(2,3)*dgmetds(1,1)+gmet(1,2)&
&       *(3*gmet(1,3)*dgmetds(1,1)+6*gmet(1,1)*dgmetds(1,3)))
       cm(6,6,3,3)=gmet(1,1)*(3*gmet(1,2)**2*dgmetds(1,1)+6*gmet(1,1)&
&       *gmet(2,2)*dgmetds(1,1)+6*gmet(1,1)*gmet(1,2)*dgmetds(1,2))
       cm(7,6,3,3)=3*gmet(1,2)**3*dgmetds(1,2)+15*gmet(1,1)*gmet(1,2)&
&       *gmet(2,2)*dgmetds(1,2)+gmet(1,2)**2*(4.5d0*gmet(2,2)*dgmetds(1,1)&
&       +3*gmet(1,1)*dgmetds(2,2))+gmet(1,1)*gmet(2,2)*(-1.5d0*gmet(2,2)&
&       *dgmetds(1,1)+6*gmet(1,1)*dgmetds(2,2))
       cm(8,6,3,3)=-3*gmet(1,2)**2*gmet(3,3)*dgmetds(1,1)+gmet(1,3)**2*(7.5d0*gmet(2,2)&
&       *dgmetds(1,1)+3*gmet(1,2)*dgmetds(1,2))+gmet(1,3)*(6*gmet(1,2)&
&       **2*dgmetds(1,3)+24*gmet(1,1)*(gmet(2,3)*dgmetds(1,2)+gmet(2,2)&
&       *dgmetds(1,3))+gmet(1,2)*(9*gmet(2,3)*dgmetds(1,1)+6*gmet(1,1)&
&       *dgmetds(2,3)))+gmet(1,1)*(-3*gmet(2,3)**2*dgmetds(1,1)-1.5d0*gmet(2,2)&
&       *gmet(3,3)*dgmetds(1,1)-9*gmet(1,2)*gmet(3,3)*dgmetds(1,2)+6*gmet(1,2)&
&       *gmet(2,3)*dgmetds(1,3)+3*gmet(1,2)**2*dgmetds(3,3))+gmet(1,1)&
&       **2*(12*gmet(2,3)*dgmetds(2,3)+6*gmet(2,2)*dgmetds(3,3))
       cm(9,6,3,3)=3*gmet(1,2)**3*dgmetds(1,3)+gmet(1,2)*(gmet(1,1)*(6*gmet(2,3)&
&       *dgmetds(1,2)+15*gmet(2,2)*dgmetds(1,3))+gmet(1,3)*(12*gmet(2,2)&
&       *dgmetds(1,1)+3*gmet(1,1)*dgmetds(2,2)))+gmet(1,2)**2*(1.5d0*gmet(2,3)&
&       *dgmetds(1,1)+6*(gmet(1,3)*dgmetds(1,2)+gmet(1,1)*dgmetds(2,3)))&
&       +gmet(1,1)*(6*gmet(1,1)*gmet(2,3)*dgmetds(2,2)+gmet(2,2)*(-4.5d0*gmet(2,3)&
&       *dgmetds(1,1)+24*gmet(1,3)*dgmetds(1,2)+12*gmet(1,1)*dgmetds(2,3)))
       cm(10,6,3,3)=gmet(1,3)**2*(7.5d0*gmet(2,3)*dgmetds(1,1)+3*gmet(1,2)&
&       *dgmetds(1,3))+gmet(1,3)*(24*gmet(1,1)*gmet(2,3)*dgmetds(1,3)&
&       +gmet(1,2)*(-3*gmet(3,3)*dgmetds(1,1)+3*gmet(1,1)*dgmetds(3,3)))&
&       +gmet(1,1)*(-9*gmet(1,2)*gmet(3,3)*dgmetds(1,3)+gmet(2,3)*(-1.5d0*gmet(3,3)&
&       *dgmetds(1,1)+6*gmet(1,1)*dgmetds(3,3)))
       cm(11,6,3,3)=9*gmet(1,2)**2*gmet(2,2)*dgmetds(1,2)-3*gmet(1,1)&
&       *gmet(2,2)**2*dgmetds(1,2)+1.5d0*gmet(1,2)**3*dgmetds(2,2)+7.5d0*gmet(1,1)&
&       *gmet(1,2)*gmet(2,2)*dgmetds(2,2)
       cm(12,6,3,3)=gmet(1,3)**2*(15*gmet(2,2)*dgmetds(1,2)+1.5d0*gmet(1,2)&
&       *dgmetds(2,2))+gmet(1,3)*(gmet(1,2)*(18*gmet(2,3)*dgmetds(1,2)&
&       +24*gmet(2,2)*dgmetds(1,3))+6*gmet(1,2)**2*dgmetds(2,3)+gmet(1,1)&
&       *(12*gmet(2,3)*dgmetds(2,2)+24*gmet(2,2)*dgmetds(2,3)))+gmet(1,2)&
&       **2*(-6*gmet(3,3)*dgmetds(1,2)+3*gmet(2,3)*dgmetds(1,3)+1.5d0*gmet(1,2)&
&       *dgmetds(3,3))+gmet(1,1)*(-6*gmet(2,3)**2*dgmetds(1,2)-3*gmet(2,2)&
&       *gmet(3,3)*dgmetds(1,2)-9*gmet(2,2)*gmet(2,3)*dgmetds(1,3)-4.5d0*gmet(1,2)&
&       *gmet(3,3)*dgmetds(2,2)+6*gmet(1,2)*gmet(2,3)*dgmetds(2,3)+7.5d0*gmet(1,2)&
&       *gmet(2,2)*dgmetds(3,3))
       cm(13,6,3,3)=gmet(1,2)**2*(3*gmet(2,3)*dgmetds(1,2)+9*gmet(2,2)&
&       *dgmetds(1,3)+3*gmet(1,3)*dgmetds(2,2))+gmet(1,1)*gmet(2,2)*(-9*gmet(2,3)&
&       *dgmetds(1,2)-3*gmet(2,2)*dgmetds(1,3)+12*gmet(1,3)*dgmetds(2,2))&
&       +3*gmet(1,2)**3*dgmetds(2,3)+gmet(1,2)*(24*gmet(1,3)*gmet(2,2)&
&       *dgmetds(1,2)+gmet(1,1)*(3*gmet(2,3)*dgmetds(2,2)+15*gmet(2,2)&
&       *dgmetds(2,3)))
       cm(14,6,3,3)=-6*gmet(1,2)**2*gmet(3,3)*dgmetds(1,3)+gmet(1,3)&
&       **2*(15*gmet(2,3)*dgmetds(1,2)+15*gmet(2,2)*dgmetds(1,3)+3*gmet(1,2)&
&       *dgmetds(2,3))+gmet(1,1)*(-3*gmet(2,3)*gmet(3,3)*dgmetds(1,2)&
&       -6*gmet(2,3)**2*dgmetds(1,3)-3*gmet(2,2)*gmet(3,3)*dgmetds(1,3)&
&       -9*gmet(1,2)*gmet(3,3)*dgmetds(2,3)+3*gmet(1,2)*gmet(2,3)*dgmetds(3,3))&
&       +gmet(1,3)*(gmet(1,2)*(-6*gmet(3,3)*dgmetds(1,2)+18*gmet(2,3)&
&       *dgmetds(1,3))+3*gmet(1,2)**2*dgmetds(3,3)+gmet(1,1)*(24*gmet(2,3)&
&       *dgmetds(2,3)+12*gmet(2,2)*dgmetds(3,3)))
       cm(15,6,3,3)=gmet(1,1)*gmet(3,3)*(-3*gmet(2,3)*dgmetds(1,3)-4.5d0*gmet(1,2)&
&       *dgmetds(3,3))+gmet(1,3)**2*(15*gmet(2,3)*dgmetds(1,3)+1.5d0*gmet(1,2)&
&       *dgmetds(3,3))+gmet(1,3)*(-6*gmet(1,2)*gmet(3,3)*dgmetds(1,3)&
&       +12*gmet(1,1)*gmet(2,3)*dgmetds(3,3))
       cm(16,6,3,3)=gmet(2,2)*(4.5d0*gmet(1,2)**2-1.5d0*gmet(1,1)*gmet(2,2))&
&       *dgmetds(2,2)
       cm(17,6,3,3)=7.5d0*gmet(1,3)**2*gmet(2,2)*dgmetds(2,2)+gmet(1,2)&
&       *gmet(1,3)*(9*gmet(2,3)*dgmetds(2,2)+24*gmet(2,2)*dgmetds(2,3))&
&       +gmet(1,2)**2*(-3*gmet(3,3)*dgmetds(2,2)+3*gmet(2,3)*dgmetds(2,3)&
&       +4.5d0*gmet(2,2)*dgmetds(3,3))+gmet(1,1)*(-3*gmet(2,3)**2*dgmetds(2,2)&
&       -9*gmet(2,2)*gmet(2,3)*dgmetds(2,3)-1.5d0*gmet(2,2)*(gmet(3,3)&
&       *dgmetds(2,2)+gmet(2,2)*dgmetds(3,3)))
       cm(18,6,3,3)=12*gmet(1,2)*gmet(1,3)*gmet(2,2)*dgmetds(2,2)+gmet(1,1)&
&       *gmet(2,2)*(-4.5d0*gmet(2,3)*dgmetds(2,2)-3*gmet(2,2)*dgmetds(2,3))&
&       +gmet(1,2)**2*(1.5d0*gmet(2,3)*dgmetds(2,2)+9*gmet(2,2)*dgmetds(2,3))
       cm(19,6,3,3)=gmet(1,3)**2*(7.5d0*gmet(2,3)*dgmetds(2,2)+15*gmet(2,2)&
&       *dgmetds(2,3))+gmet(1,2)*gmet(1,3)*(-3*gmet(3,3)*dgmetds(2,2)&
&       +18*gmet(2,3)*dgmetds(2,3)+12*gmet(2,2)*dgmetds(3,3))+gmet(1,2)&
&       **2*(-6*gmet(3,3)*dgmetds(2,3)+1.5d0*gmet(2,3)*dgmetds(3,3))&
&       +gmet(1,1)*(-1.5d0*gmet(2,3)*gmet(3,3)*dgmetds(2,2)-6*gmet(2,3)&
&       **2*dgmetds(2,3)-3*gmet(2,2)*gmet(3,3)*dgmetds(2,3)-4.5d0*gmet(2,2)&
&       *gmet(2,3)*dgmetds(3,3))
       cm(20,6,3,3)=(8*(90*gmet(1,3)**2*gmet(2,3)-36*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)-18*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(2,3)+12*(30*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(-12*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(3,3))&
&       /48.d0
       cm(21,6,3,3)=((90*gmet(1,3)**2*gmet(2,3)-36*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)-18*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(3,3))/12.d0
       cm(1,7,3,3)=gmet(1,2)*(2.5d0*gmet(1,2)**2-1.5d0*gmet(1,1)*gmet(2,2))&
&       *dgmetds(1,1)
       cm(2,7,3,3)=9*gmet(1,2)**2*gmet(2,2)*dgmetds(1,2)-3*gmet(1,1)&
&       *gmet(2,2)**2*dgmetds(1,2)+2.5d0*gmet(1,2)**3*dgmetds(2,2)+gmet(1,2)&
&       *gmet(2,2)*(3*gmet(2,2)*dgmetds(1,1)-1.5d0*gmet(1,1)*dgmetds(2,2))
       cm(3,7,3,3)=(6*(-36*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(90*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+12*(-36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,2)+90*gmet(1,2)**2*gmet(2,3)-18*gmet(1,1)*gmet(2,2)*gmet(2,3))&
&       *dgmetds(1,3)+2*(90*gmet(1,2)**3-54*gmet(1,1)*gmet(1,2)*gmet(2,2))&
&       *dgmetds(3,3))/72.d0
       cm(4,7,3,3)=(6*gmet(2,2)*(-18*gmet(1,3)*gmet(2,2)+54*gmet(1,2)&
&       *gmet(2,3))*dgmetds(1,1)+6*(-36*gmet(1,2)*gmet(1,3)*gmet(2,2)&
&       +90*gmet(1,2)**2*gmet(2,3)-18*gmet(1,1)*gmet(2,2)*gmet(2,3))&
&       *dgmetds(1,2)+6*gmet(2,2)*(54*gmet(1,2)**2-18*gmet(1,1)*gmet(2,2))&
&       *dgmetds(1,3)+2*(90*gmet(1,2)**3-54*gmet(1,1)*gmet(1,2)*gmet(2,2))&
&       *dgmetds(2,3))/36.d0
       cm(5,7,3,3)=7.5d0*gmet(1,2)**2*gmet(2,3)*dgmetds(1,1)-1.5d0*gmet(1,1)&
&       *gmet(2,2)*gmet(2,3)*dgmetds(1,1)+5*gmet(1,2)**3*dgmetds(1,3)&
&       -3*gmet(1,2)*gmet(2,2)*(gmet(1,3)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,3))
       cm(6,7,3,3)=4.5d0*gmet(1,2)**2*gmet(2,2)*dgmetds(1,1)-1.5d0*gmet(1,1)&
&       *gmet(2,2)**2*dgmetds(1,1)+5*gmet(1,2)**3*dgmetds(1,2)-3*gmet(1,1)&
&       *gmet(1,2)*gmet(2,2)*dgmetds(1,2)
       cm(7,7,3,3)=gmet(2,2)*(1*gmet(2,2)**2*dgmetds(1,1)+4.5d0*gmet(1,2)&
&       **2*dgmetds(2,2)+gmet(2,2)*(6*gmet(1,2)*dgmetds(1,2)-1.5d0*gmet(1,1)&
&       *dgmetds(2,2)))
       cm(8,7,3,3)=(180*gmet(1,2)*gmet(2,3)*(gmet(2,3)*dgmetds(1,2)+gmet(1,2)&
&       *dgmetds(2,3))+gmet(2,2)**2*(-18*gmet(3,3)*dgmetds(1,1)-72*gmet(1,3)&
&       *dgmetds(1,3)-18*gmet(1,1)*dgmetds(3,3))+gmet(2,2)*(54*gmet(2,3)&
&       **2*dgmetds(1,1)+gmet(2,3)*(-72*gmet(1,3)*dgmetds(1,2)+216*gmet(1,2)&
&       *dgmetds(1,3)-36*gmet(1,1)*dgmetds(2,3))+gmet(1,2)*(-36*gmet(3,3)&
&       *dgmetds(1,2)-72*gmet(1,3)*dgmetds(2,3)+54*gmet(1,2)*dgmetds(3,3))))&
&       /12.d0
       cm(9,7,3,3)=(180*gmet(1,2)**2*gmet(2,3)*dgmetds(2,2)+gmet(2,2)&
&       **2*(72*gmet(2,3)*dgmetds(1,1)-144*gmet(1,3)*dgmetds(1,2)+144*gmet(1,2)&
&       *dgmetds(1,3)-72*gmet(1,1)*dgmetds(2,3))+gmet(2,2)*(-36*gmet(1,1)&
&       *gmet(2,3)*dgmetds(2,2)+gmet(1,2)*(432*gmet(2,3)*dgmetds(1,2)&
&       -72*gmet(1,3)*dgmetds(2,2))+216*gmet(1,2)**2*dgmetds(2,3)))/24.d0
       cm(10,7,3,3)=(180*gmet(2,3)**3*dgmetds(1,1)+1080*gmet(1,2)*gmet(2,3)&
&       **2*dgmetds(1,3)-216*gmet(1,2)*gmet(2,2)*(gmet(3,3)*dgmetds(1,3)&
&       +gmet(1,3)*dgmetds(3,3))+gmet(2,3)*(540*gmet(1,2)**2*dgmetds(3,3)&
&       +gmet(2,2)*(-108*gmet(3,3)*dgmetds(1,1)-432*gmet(1,3)*dgmetds(1,3)&
&       -108*gmet(1,1)*dgmetds(3,3))))/72.d0
       cm(11,7,3,3)=gmet(2,2)**2*(2*gmet(2,2)*dgmetds(1,2)+3*gmet(1,2)&
&       *dgmetds(2,2))
       cm(12,7,3,3)=(180*gmet(1,2)*gmet(2,3)**2*dgmetds(2,2)+gmet(2,2)&
&       *(216*gmet(2,3)**2*dgmetds(1,2)-36*gmet(1,2)*gmet(3,3)*dgmetds(2,2)&
&       +gmet(2,3)*(-72*gmet(1,3)*dgmetds(2,2)+432*gmet(1,2)*dgmetds(2,3)))&
&       +gmet(2,2)**2*(-72*gmet(3,3)*dgmetds(1,2)+144*gmet(2,3)*dgmetds(1,3)&
&       -144*gmet(1,3)*dgmetds(2,3)+72*gmet(1,2)*dgmetds(3,3)))/24.d0
       cm(13,7,3,3)=gmet(2,2)*(2*gmet(2,2)**2*dgmetds(1,3)+9*gmet(1,2)&
&       *gmet(2,3)*dgmetds(2,2)+gmet(2,2)*(6*gmet(2,3)*dgmetds(1,2)-3*gmet(1,3)&
&       *dgmetds(2,2)+6*gmet(1,2)*dgmetds(2,3)))
       cm(14,7,3,3)=(180*gmet(2,3)**3*dgmetds(1,2)+gmet(2,3)**2*(324*gmet(2,2)&
&       *dgmetds(1,3)+540*gmet(1,2)*dgmetds(2,3))+gmet(2,2)*gmet(2,3)&
&       *(-108*gmet(3,3)*dgmetds(1,2)-216*gmet(1,3)*dgmetds(2,3)+324*gmet(1,2)&
&       *dgmetds(3,3))-108*gmet(2,2)*(gmet(1,2)*gmet(3,3)*dgmetds(2,3)&
&       +gmet(2,2)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))))&
&       /36.d0
       cm(15,7,3,3)=5*gmet(2,3)**3*dgmetds(1,3)+7.5d0*gmet(1,2)*gmet(2,3)&
&       **2*dgmetds(3,3)-1.5d0*gmet(1,2)*gmet(2,2)*gmet(3,3)*dgmetds(3,3)&
&       -3*gmet(2,2)*gmet(2,3)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)*dgmetds(3,3))
       cm(16,7,3,3)=gmet(2,2)**3*dgmetds(2,2)
       cm(17,7,3,3)=gmet(2,2)*(4.5d0*gmet(2,3)**2*dgmetds(2,2)+6*gmet(2,2)&
&       *gmet(2,3)*dgmetds(2,3)+gmet(2,2)*(-1.5d0*gmet(3,3)*dgmetds(2,2)&
&       +gmet(2,2)*dgmetds(3,3)))
       cm(18,7,3,3)=gmet(2,2)**2*(3*gmet(2,3)*dgmetds(2,2)+2*gmet(2,2)&
&       *dgmetds(2,3))
       cm(19,7,3,3)=2.5d0*gmet(2,3)**3*dgmetds(2,2)+9*gmet(2,2)*gmet(2,3)&
&       **2*dgmetds(2,3)-3*gmet(2,2)**2*gmet(3,3)*dgmetds(2,3)+gmet(2,2)&
&       *gmet(2,3)*(-1.5d0*gmet(3,3)*dgmetds(2,2)+3*gmet(2,2)*dgmetds(3,3))
       cm(20,7,3,3)=5*gmet(2,3)**3*dgmetds(2,3)-3*gmet(2,2)*gmet(2,3)&
&       *gmet(3,3)*dgmetds(2,3)+4.5d0*gmet(2,2)*gmet(2,3)**2*dgmetds(3,3)&
&       -1.5d0*gmet(2,2)**2*gmet(3,3)*dgmetds(3,3)
       cm(21,7,3,3)=gmet(2,3)*(2.5d0*gmet(2,3)**2-1.5d0*gmet(2,2)*gmet(3,3))&
&       *dgmetds(3,3)
       cm(1,8,3,3)=((-36*gmet(1,1)*gmet(1,3)*gmet(2,3)+gmet(1,2)*(90*gmet(1,3)&
&       **2-18*gmet(1,1)*gmet(3,3)))*dgmetds(1,1))/12.d0
       cm(2,8,3,3)=(6*(48*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)*(6*gmet(2,3)&
&       **2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+12*(30*gmet(1,3)**2*gmet(2,2)&
&       +36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)&
&       *(-12*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+2*(-36*gmet(1,1)&
&       *gmet(1,3)*gmet(2,3)+gmet(1,2)*(90*gmet(1,3)**2-18*gmet(1,1)&
&       *gmet(3,3)))*dgmetds(2,2))/24.d0
       cm(3,8,3,3)=(6*gmet(3,3)*(12*gmet(1,3)*gmet(2,3)+24*gmet(1,2)&
&       *gmet(3,3))*dgmetds(1,1)+12*(6*gmet(1,3)**2*gmet(2,3)+48*gmet(1,2)&
&       *gmet(1,3)*gmet(3,3)-18*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(1,3)&
&       +2*(-36*gmet(1,1)*gmet(1,3)*gmet(2,3)+gmet(1,2)*(90*gmet(1,3)&
&       **2-18*gmet(1,1)*gmet(3,3)))*dgmetds(3,3))/24.d0
       cm(4,8,3,3)=(6*(6*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(6*gmet(2,3)&
&       **2+24*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+6*(6*gmet(1,3)**2*gmet(2,3)&
&       +48*gmet(1,2)*gmet(1,3)*gmet(3,3)-18*gmet(1,1)*gmet(2,3)*gmet(3,3))&
&       *dgmetds(1,2)+6*(30*gmet(1,3)**2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(-12*gmet(2,3)&
&       **2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,3)+2*(-36*gmet(1,1)*gmet(1,3)&
&       *gmet(2,3)+gmet(1,2)*(90*gmet(1,3)**2-18*gmet(1,1)*gmet(3,3)))&
&       *dgmetds(2,3))/12.d0
       cm(5,8,3,3)=(12*(6*gmet(1,3)**2*gmet(2,3)+48*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)-18*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(1,1)+8*(-36*gmet(1,1)&
&       *gmet(1,3)*gmet(2,3)+gmet(1,2)*(90*gmet(1,3)**2-18*gmet(1,1)&
&       *gmet(3,3)))*dgmetds(1,3))/48.d0
       cm(6,8,3,3)=(12*(30*gmet(1,3)**2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(-12*gmet(2,3)&
&       **2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+8*(-36*gmet(1,1)*gmet(1,3)&
&       *gmet(2,3)+gmet(1,2)*(90*gmet(1,3)**2-18*gmet(1,1)*gmet(3,3)))&
&       *dgmetds(1,2))/48.d0
       cm(7,8,3,3)=(2*gmet(2,2)*(54*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3))&
&       *dgmetds(1,1)+12*(48*gmet(1,3)*gmet(2,2)*gmet(2,3)+gmet(1,2)&
&       *(6*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+6*(30*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(-12*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(2,2))&
&       /24.d0
       cm(8,8,3,3)=(gmet(3,3)*(12*gmet(2,3)**2+24*gmet(2,2)*gmet(3,3))&
&       *dgmetds(1,1)+2*gmet(3,3)*(12*gmet(1,3)*gmet(2,3)+24*gmet(1,2)&
&       *gmet(3,3))*dgmetds(1,2)+4*(6*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)&
&       *(6*gmet(2,3)**2+24*gmet(2,2)*gmet(3,3)))*dgmetds(1,3)+2*(6*gmet(1,3)&
&       **2*gmet(2,3)+48*gmet(1,2)*gmet(1,3)*gmet(3,3)-18*gmet(1,1)*gmet(2,3)&
&       *gmet(3,3))*dgmetds(2,3)+(30*gmet(1,3)**2*gmet(2,2)+36*gmet(1,2)&
&       *gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(-12*gmet(2,3)&
&       **2-6*gmet(2,2)*gmet(3,3)))*dgmetds(3,3))/4.d0
       cm(9,8,3,3)=((6*gmet(2,3)**3+30*gmet(2,2)*gmet(2,3)*gmet(3,3))&
&       *dgmetds(1,1)+4*(6*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)*(6*gmet(2,3)&
&       **2+24*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+2*(48*gmet(1,3)*gmet(2,2)&
&       *gmet(2,3)+gmet(1,2)*(6*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))&
&       *dgmetds(1,3)+(6*gmet(1,3)**2*gmet(2,3)+48*gmet(1,2)*gmet(1,3)&
&       *gmet(3,3)-18*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(2,2)+2*(30*gmet(1,3)&
&       **2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)*gmet(2,3)-12*gmet(1,2)**2*gmet(3,3)&
&       +gmet(1,1)*(-12*gmet(2,3)**2-6*gmet(2,2)*gmet(3,3)))*dgmetds(2,3))&
&       /4.d0
       cm(10,8,3,3)=(288*gmet(1,2)*gmet(3,3)*(gmet(3,3)*dgmetds(1,3)&
&       +gmet(1,3)*dgmetds(3,3))+gmet(2,3)*(72*gmet(3,3)**2*dgmetds(1,1)&
&       +36*gmet(1,3)**2*dgmetds(3,3)+gmet(3,3)*(144*gmet(1,3)*dgmetds(1,3)&
&       -108*gmet(1,1)*dgmetds(3,3))))/24.d0
       cm(11,8,3,3)=-3*gmet(2,2)**2*gmet(3,3)*dgmetds(1,2)+1.5d0*gmet(1,2)&
&       *gmet(2,3)**2*dgmetds(2,2)+gmet(2,2)*(9*gmet(2,3)**2*dgmetds(1,2)&
&       +12*gmet(1,3)*gmet(2,3)*dgmetds(2,2)-4.5d0*gmet(1,2)*gmet(3,3)&
&       *dgmetds(2,2))
       cm(12,8,3,3)=(2*gmet(3,3)*(12*gmet(2,3)**2+24*gmet(2,2)*gmet(3,3))&
&       *dgmetds(1,2)+2*(6*gmet(2,3)**3+30*gmet(2,2)*gmet(2,3)*gmet(3,3))&
&       *dgmetds(1,3)+gmet(3,3)*(12*gmet(1,3)*gmet(2,3)+24*gmet(1,2)&
&       *gmet(3,3))*dgmetds(2,2)+4*(6*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)&
&       *(6*gmet(2,3)**2+24*gmet(2,2)*gmet(3,3)))*dgmetds(2,3)+(48*gmet(1,3)&
&       *gmet(2,2)*gmet(2,3)+gmet(1,2)*(6*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))&
&       *dgmetds(3,3))/4.d0
       cm(13,8,3,3)=(36*gmet(2,3)**3*dgmetds(1,2)+gmet(2,2)*gmet(3,3)&
&       *(-36*gmet(2,2)*dgmetds(1,3)+144*gmet(1,3)*dgmetds(2,2)-108*gmet(1,2)&
&       *dgmetds(2,3))+gmet(2,3)**2*(108*gmet(2,2)*dgmetds(1,3)+36*(gmet(1,3)&
&       *dgmetds(2,2)+gmet(1,2)*dgmetds(2,3)))+gmet(2,3)*(36*gmet(1,2)&
&       *gmet(3,3)*dgmetds(2,2)+gmet(2,2)*(180*gmet(3,3)*dgmetds(1,2)&
&       +288*gmet(1,3)*dgmetds(2,3))))/12.d0
       cm(14,8,3,3)=(gmet(2,3)*gmet(3,3)*(72*gmet(3,3)*dgmetds(1,2)+72*gmet(1,3)&
&       *dgmetds(2,3)+36*gmet(1,2)*dgmetds(3,3))+gmet(2,3)**2*(72*gmet(3,3)&
&       *dgmetds(1,3)+36*gmet(1,3)*dgmetds(3,3))+144*gmet(3,3)*(gmet(1,2)&
&       *gmet(3,3)*dgmetds(2,3)+gmet(2,2)*(gmet(3,3)*dgmetds(1,3)+gmet(1,3)&
&       *dgmetds(3,3))))/12.d0
       cm(15,8,3,3)=gmet(3,3)*(6*gmet(1,2)*gmet(3,3)*dgmetds(3,3)+gmet(2,3)&
&       *(6*gmet(3,3)*dgmetds(1,3)+3*gmet(1,3)*dgmetds(3,3)))
       cm(16,8,3,3)=gmet(2,2)*(4.5d0*gmet(2,3)**2-1.5d0*gmet(2,2)*gmet(3,3))&
&       *dgmetds(2,2)
       cm(17,8,3,3)=3*gmet(2,3)**3*dgmetds(2,3)+15*gmet(2,2)*gmet(2,3)&
&       *gmet(3,3)*dgmetds(2,3)+gmet(2,2)*gmet(3,3)*(6*gmet(3,3)*dgmetds(2,2)&
&       -1.5d0*gmet(2,2)*dgmetds(3,3))+gmet(2,3)**2*(3*gmet(3,3)*dgmetds(2,2)&
&       +4.5d0*gmet(2,2)*dgmetds(3,3))
       cm(18,8,3,3)=1.5d0*gmet(2,3)**3*dgmetds(2,2)+7.5d0*gmet(2,2)*gmet(2,3)&
&       *gmet(3,3)*dgmetds(2,2)+9*gmet(2,2)*gmet(2,3)**2*dgmetds(2,3)&
&       -3*gmet(2,2)**2*gmet(3,3)*dgmetds(2,3)
       cm(19,8,3,3)=6*gmet(2,3)**2*gmet(3,3)*dgmetds(2,3)+12*gmet(2,2)&
&       *gmet(3,3)**2*dgmetds(2,3)+1.5d0*gmet(2,3)**3*dgmetds(3,3)+gmet(2,3)&
&       *gmet(3,3)*(3*gmet(3,3)*dgmetds(2,2)+7.5d0*gmet(2,2)*dgmetds(3,3))
       cm(20,8,3,3)=gmet(3,3)*(6*gmet(2,3)*gmet(3,3)*dgmetds(2,3)+3*gmet(2,3)&
&       **2*dgmetds(3,3)+6*gmet(2,2)*gmet(3,3)*dgmetds(3,3))
       cm(21,8,3,3)=3*gmet(2,3)*gmet(3,3)**2*dgmetds(3,3)
       cm(1,9,3,3)=((90*gmet(1,2)**2*gmet(1,3)-18*gmet(1,1)*gmet(1,3)&
&       *gmet(2,2)-36*gmet(1,1)*gmet(1,2)*gmet(2,3))*dgmetds(1,1))/12.d0
       cm(2,9,3,3)=gmet(2,3)*(3*gmet(1,2)**2*dgmetds(1,2)-9*gmet(1,1)&
&       *gmet(2,2)*dgmetds(1,2)+gmet(1,2)*(3*gmet(2,2)*dgmetds(1,1)-3*gmet(1,1)&
&       *dgmetds(2,2)))+gmet(1,3)*(6*gmet(2,2)**2*dgmetds(1,1)+7.5d0*gmet(1,2)&
&       **2*dgmetds(2,2)+gmet(2,2)*(24*gmet(1,2)*dgmetds(1,2)-1.5d0*gmet(1,1)&
&       *dgmetds(2,2)))
       cm(3,9,3,3)=-6*gmet(1,3)**2*gmet(2,2)*dgmetds(1,3)+15*gmet(1,2)&
&       **2*gmet(3,3)*dgmetds(1,3)+gmet(1,1)*(-6*gmet(2,3)**2-3*gmet(2,2)&
&       *gmet(3,3))*dgmetds(1,3)+gmet(1,2)*gmet(2,3)*(12*gmet(3,3)*dgmetds(1,1)&
&       -3*gmet(1,1)*dgmetds(3,3))+gmet(1,3)*(1.5d0*gmet(2,3)**2*dgmetds(1,1)&
&       -4.5d0*gmet(2,2)*gmet(3,3)*dgmetds(1,1)+18*gmet(1,2)*gmet(2,3)&
&       *dgmetds(1,3)+7.5d0*gmet(1,2)**2*dgmetds(3,3)-1.5d0*gmet(1,1)&
&       *gmet(2,2)*dgmetds(3,3))
       cm(4,9,3,3)=-6*gmet(1,3)**2*gmet(2,2)*dgmetds(1,2)+gmet(1,2)**2*(15*gmet(3,3)&
&       *dgmetds(1,2)+3*gmet(2,3)*dgmetds(1,3))+gmet(1,1)*(-6*gmet(2,3)&
&       **2*dgmetds(1,2)-3*gmet(2,2)*gmet(3,3)*dgmetds(1,2)-9*gmet(2,2)&
&       *gmet(2,3)*dgmetds(1,3))+gmet(1,2)*(3*gmet(2,3)**2*dgmetds(1,1)&
&       +12*gmet(2,2)*gmet(3,3)*dgmetds(1,1)-6*gmet(1,1)*gmet(2,3)*dgmetds(2,3))&
&       +gmet(1,3)*(gmet(2,2)*(3*gmet(2,3)*dgmetds(1,1)+24*gmet(1,2)&
&       *dgmetds(1,3)-3*gmet(1,1)*dgmetds(2,3))+gmet(1,2)*(18*gmet(2,3)&
&       *dgmetds(1,2)+15*gmet(1,2)*dgmetds(2,3)))
       cm(5,9,3,3)=(12*(-12*gmet(1,3)**2*gmet(2,2)+36*gmet(1,2)*gmet(1,3)&
&       *gmet(2,3)+30*gmet(1,2)**2*gmet(3,3)+gmet(1,1)*(-12*gmet(2,3)&
&       **2-6*gmet(2,2)*gmet(3,3)))*dgmetds(1,1)+8*(90*gmet(1,2)**2*gmet(1,3)&
&       -18*gmet(1,1)*gmet(1,3)*gmet(2,2)-36*gmet(1,1)*gmet(1,2)*gmet(2,3))&
&       *dgmetds(1,3))/48.d0
       cm(6,9,3,3)=gmet(1,1)*gmet(2,2)*(-4.5d0*gmet(2,3)*dgmetds(1,1)&
&       -3*gmet(1,3)*dgmetds(1,2))+gmet(1,2)**2*(1.5d0*gmet(2,3)*dgmetds(1,1)&
&       +15*gmet(1,3)*dgmetds(1,2))+gmet(1,2)*(12*gmet(1,3)*gmet(2,2)&
&       *dgmetds(1,1)-6*gmet(1,1)*gmet(2,3)*dgmetds(1,2))
       cm(7,9,3,3)=gmet(2,2)**2*(3*gmet(2,3)*dgmetds(1,1)+12*gmet(1,3)&
&       *dgmetds(1,2))+1.5d0*gmet(1,2)**2*gmet(2,3)*dgmetds(2,2)+gmet(2,2)&
&       *(-4.5d0*gmet(1,1)*gmet(2,3)*dgmetds(2,2)+gmet(1,2)*(6*gmet(2,3)&
&       *dgmetds(1,2)+12*gmet(1,3)*dgmetds(2,2)))
       cm(8,9,3,3)=1.5d0*gmet(2,3)**3*dgmetds(1,1)-6*gmet(1,3)**2*gmet(2,2)&
&       *dgmetds(2,3)+gmet(2,3)**2*(3*gmet(1,3)*dgmetds(1,2)+6*gmet(1,2)&
&       *dgmetds(1,3)-6*gmet(1,1)*dgmetds(2,3))+gmet(3,3)*(24*gmet(1,2)&
&       *gmet(2,2)*dgmetds(1,3)+15*gmet(1,2)**2*dgmetds(2,3)-3*gmet(1,1)&
&       *gmet(2,2)*dgmetds(2,3))+gmet(1,3)*gmet(2,2)*(-9*gmet(3,3)*dgmetds(1,2)&
&       +12*gmet(1,2)*dgmetds(3,3))+gmet(2,3)*(gmet(2,2)*(7.5d0*gmet(3,3)&
&       *dgmetds(1,1)+6*gmet(1,3)*dgmetds(1,3)-4.5d0*gmet(1,1)*dgmetds(3,3))&
&       +gmet(1,2)*(24*gmet(3,3)*dgmetds(1,2)+18*gmet(1,3)*dgmetds(2,3)&
&       +1.5d0*gmet(1,2)*dgmetds(3,3)))
       cm(9,9,3,3)=gmet(2,2)**2*(6*gmet(3,3)*dgmetds(1,1)+12*gmet(1,3)&
&       *dgmetds(1,3))-3*gmet(1,1)*gmet(2,3)**2*dgmetds(2,2)+gmet(1,2)&
&       *gmet(2,3)*(6*gmet(2,3)*dgmetds(1,2)+9*gmet(1,3)*dgmetds(2,2))&
&       +gmet(1,2)**2*(7.5d0*gmet(3,3)*dgmetds(2,2)+3*gmet(2,3)*dgmetds(2,3))&
&       +gmet(2,2)*(3*gmet(2,3)**2*dgmetds(1,1)+(-3*gmet(1,3)**2-1.5d0*gmet(1,1)&
&       *gmet(3,3))*dgmetds(2,2)+gmet(2,3)*(6*gmet(1,3)*dgmetds(1,2)&
&       +6*gmet(1,2)*dgmetds(1,3)-9*gmet(1,1)*dgmetds(2,3))+24*gmet(1,2)&
&       *(gmet(3,3)*dgmetds(1,2)+gmet(1,3)*dgmetds(2,3)))
       cm(10,9,3,3)=7.5d0*gmet(1,2)**2*gmet(3,3)*dgmetds(3,3)+gmet(2,3)&
&       **2*(4.5d0*gmet(3,3)*dgmetds(1,1)+3*gmet(1,3)*dgmetds(1,3)-3*gmet(1,1)&
&       *dgmetds(3,3))+gmet(1,2)*gmet(2,3)*(24*gmet(3,3)*dgmetds(1,3)&
&       +9*gmet(1,3)*dgmetds(3,3))+gmet(2,2)*(-1.5d0*gmet(3,3)**2*dgmetds(1,1)&
&       -9*gmet(1,3)*gmet(3,3)*dgmetds(1,3)-3*gmet(1,3)**2*dgmetds(3,3)&
&       -1.5d0*gmet(1,1)*gmet(3,3)*dgmetds(3,3))
       cm(11,9,3,3)=gmet(2,2)*(3*gmet(1,2)*gmet(2,3)*dgmetds(2,2)+6*gmet(2,2)&
&       *(gmet(2,3)*dgmetds(1,2)+gmet(1,3)*dgmetds(2,2)))
       cm(12,9,3,3)=3*gmet(2,3)**3*dgmetds(1,2)+gmet(2,3)**2*(6*gmet(2,2)&
&       *dgmetds(1,3)+1.5d0*gmet(1,3)*dgmetds(2,2)+6*gmet(1,2)*dgmetds(2,3))&
&       +gmet(2,3)*(12*gmet(1,2)*gmet(3,3)*dgmetds(2,2)+gmet(2,2)*(15*gmet(3,3)&
&       *dgmetds(1,2)+6*gmet(1,3)*dgmetds(2,3)+3*gmet(1,2)*dgmetds(3,3)))&
&       +gmet(2,2)*(gmet(3,3)*(-4.5d0*gmet(1,3)*dgmetds(2,2)+24*gmet(1,2)&
&       *dgmetds(2,3))+gmet(2,2)*(12*gmet(3,3)*dgmetds(1,3)+6*gmet(1,3)&
&       *dgmetds(3,3)))
       cm(13,9,3,3)=3*gmet(1,2)*gmet(2,3)**2*dgmetds(2,2)+gmet(2,2)**2*(12*gmet(3,3)&
&       *dgmetds(1,2)+6*gmet(2,3)*dgmetds(1,3)+12*gmet(1,3)*dgmetds(2,3))&
&       +gmet(2,2)*(6*gmet(2,3)**2*dgmetds(1,2)+12*gmet(1,2)*gmet(3,3)&
&       *dgmetds(2,2)+gmet(2,3)*(3*gmet(1,3)*dgmetds(2,2)+6*gmet(1,2)&
&       *dgmetds(2,3)))
       cm(14,9,3,3)=3*gmet(2,3)**3*dgmetds(1,3)+gmet(2,2)*gmet(3,3)*(-3*gmet(3,3)&
&       *dgmetds(1,2)-9*gmet(1,3)*dgmetds(2,3)+12*gmet(1,2)*dgmetds(3,3))&
&       +gmet(2,3)**2*(9*gmet(3,3)*dgmetds(1,2)+3*(gmet(1,3)*dgmetds(2,3)&
&       +gmet(1,2)*dgmetds(3,3)))+gmet(2,3)*(24*gmet(1,2)*gmet(3,3)*dgmetds(2,3)&
&       +gmet(2,2)*(15*gmet(3,3)*dgmetds(1,3)+3*gmet(1,3)*dgmetds(3,3)))
       cm(15,9,3,3)=12*gmet(1,2)*gmet(2,3)*gmet(3,3)*dgmetds(3,3)+gmet(2,2)&
&       *gmet(3,3)*(-3*gmet(3,3)*dgmetds(1,3)-4.5d0*gmet(1,3)*dgmetds(3,3))&
&       +gmet(2,3)**2*(9*gmet(3,3)*dgmetds(1,3)+1.5d0*gmet(1,3)*dgmetds(3,3))
       cm(16,9,3,3)=3*gmet(2,2)**2*gmet(2,3)*dgmetds(2,2)
       cm(17,9,3,3)=1.5d0*gmet(2,3)**3*dgmetds(2,2)+6*gmet(2,2)*gmet(2,3)&
&       **2*dgmetds(2,3)+12*gmet(2,2)**2*gmet(3,3)*dgmetds(2,3)+gmet(2,2)&
&       *gmet(2,3)*(7.5d0*gmet(3,3)*dgmetds(2,2)+3*gmet(2,2)*dgmetds(3,3))
       cm(18,9,3,3)=gmet(2,2)*(3*gmet(2,3)**2*dgmetds(2,2)+6*gmet(2,2)&
&       *gmet(3,3)*dgmetds(2,2)+6*gmet(2,2)*gmet(2,3)*dgmetds(2,3))
       cm(19,9,3,3)=3*gmet(2,3)**3*dgmetds(2,3)+15*gmet(2,2)*gmet(2,3)&
&       *gmet(3,3)*dgmetds(2,3)+gmet(2,3)**2*(4.5d0*gmet(3,3)*dgmetds(2,2)&
&       +3*gmet(2,2)*dgmetds(3,3))+gmet(2,2)*gmet(3,3)*(-1.5d0*gmet(3,3)&
&       *dgmetds(2,2)+6*gmet(2,2)*dgmetds(3,3))
       cm(20,9,3,3)=9*gmet(2,3)**2*gmet(3,3)*dgmetds(2,3)-3*gmet(2,2)&
&       *gmet(3,3)**2*dgmetds(2,3)+1.5d0*gmet(2,3)**3*dgmetds(3,3)+7.5d0*gmet(2,2)&
&       *gmet(2,3)*gmet(3,3)*dgmetds(3,3)
       cm(21,9,3,3)=gmet(3,3)*(4.5d0*gmet(2,3)**2-1.5d0*gmet(2,2)*gmet(3,3))&
&       *dgmetds(3,3)
       cm(1,10,3,3)=gmet(1,3)*(2.5d0*gmet(1,3)**2-1.5d0*gmet(1,1)*gmet(3,3))&
&       *dgmetds(1,1)
       cm(2,10,3,3)=(1080*gmet(1,3)**2*gmet(2,3)*dgmetds(1,2)-216*gmet(2,3)&
&       *gmet(3,3)*(gmet(1,2)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,2))+180*gmet(1,3)&
&       **3*dgmetds(2,2)+gmet(1,3)*(540*gmet(2,3)**2*dgmetds(1,1)+gmet(3,3)&
&       *(-108*gmet(2,2)*dgmetds(1,1)-432*gmet(1,2)*dgmetds(1,2)-108*gmet(1,1)&
&       *dgmetds(2,2))))/72.d0
       cm(3,10,3,3)=9*gmet(1,3)**2*gmet(3,3)*dgmetds(1,3)-3*gmet(1,1)&
&       *gmet(3,3)**2*dgmetds(1,3)+2.5d0*gmet(1,3)**3*dgmetds(3,3)+gmet(1,3)&
&       *gmet(3,3)*(3*gmet(3,3)*dgmetds(1,1)-1.5d0*gmet(1,1)*dgmetds(3,3))
       cm(4,10,3,3)=(gmet(1,3)**2*(324*gmet(3,3)*dgmetds(1,2)+540*gmet(2,3)&
&       *dgmetds(1,3))-108*gmet(3,3)*(gmet(1,2)*gmet(3,3)*dgmetds(1,1)&
&       +gmet(1,1)*(gmet(3,3)*dgmetds(1,2)+gmet(2,3)*dgmetds(1,3)))+180*gmet(1,3)&
&       **3*dgmetds(2,3)+gmet(1,3)*gmet(3,3)*(324*gmet(2,3)*dgmetds(1,1)&
&       -216*gmet(1,2)*dgmetds(1,3)-108*gmet(1,1)*dgmetds(2,3)))/36.d0
       cm(5,10,3,3)=4.5d0*gmet(1,3)**2*gmet(3,3)*dgmetds(1,1)-1.5d0*gmet(1,1)&
&       *gmet(3,3)**2*dgmetds(1,1)+5*gmet(1,3)**3*dgmetds(1,3)-3*gmet(1,1)&
&       *gmet(1,3)*gmet(3,3)*dgmetds(1,3)
       cm(6,10,3,3)=7.5d0*gmet(1,3)**2*gmet(2,3)*dgmetds(1,1)-1.5d0*gmet(1,1)&
&       *gmet(2,3)*gmet(3,3)*dgmetds(1,1)+5*gmet(1,3)**3*dgmetds(1,2)&
&       -3*gmet(1,3)*gmet(3,3)*(gmet(1,2)*dgmetds(1,1)+gmet(1,1)*dgmetds(1,2))
       cm(7,10,3,3)=(2*(90*gmet(2,3)**3-54*gmet(2,2)*gmet(2,3)*gmet(3,3))&
&       *dgmetds(1,1)+12*(-36*gmet(1,2)*gmet(2,3)*gmet(3,3)+gmet(1,3)&
&       *(90*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,2)+6*(90*gmet(1,3)&
&       **2*gmet(2,3)-36*gmet(1,2)*gmet(1,3)*gmet(3,3)-18*gmet(1,1)*gmet(2,3)&
&       *gmet(3,3))*dgmetds(2,2))/72.d0
       cm(8,10,3,3)=(gmet(2,3)*(72*gmet(3,3)**2*dgmetds(1,1)+180*gmet(1,3)&
&       **2*dgmetds(3,3)+gmet(3,3)*(432*gmet(1,3)*dgmetds(1,3)-36*gmet(1,1)&
&       *dgmetds(3,3)))+gmet(3,3)*(216*gmet(1,3)**2*dgmetds(2,3)+gmet(3,3)&
&       *(-144*gmet(1,2)*dgmetds(1,3)-72*gmet(1,1)*dgmetds(2,3))+gmet(1,3)&
&       *(144*gmet(3,3)*dgmetds(1,2)-72*gmet(1,2)*dgmetds(3,3))))/24.d0
       cm(9,10,3,3)=(gmet(3,3)*(54*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3))&
&       *dgmetds(1,1)+4*gmet(3,3)*(54*gmet(1,3)*gmet(2,3)-18*gmet(1,2)&
&       *gmet(3,3))*dgmetds(1,2)+2*(-36*gmet(1,2)*gmet(2,3)*gmet(3,3)&
&       +gmet(1,3)*(90*gmet(2,3)**2-18*gmet(2,2)*gmet(3,3)))*dgmetds(1,3)&
&       +gmet(3,3)*(54*gmet(1,3)**2-18*gmet(1,1)*gmet(3,3))*dgmetds(2,2)&
&       +2*(90*gmet(1,3)**2*gmet(2,3)-36*gmet(1,2)*gmet(1,3)*gmet(3,3)&
&       -18*gmet(1,1)*gmet(2,3)*gmet(3,3))*dgmetds(2,3))/12.d0
       cm(10,10,3,3)=gmet(3,3)*(1*gmet(3,3)**2*dgmetds(1,1)+4.5d0*gmet(1,3)&
&       **2*dgmetds(3,3)+gmet(3,3)*(6*gmet(1,3)*dgmetds(1,3)-1.5d0*gmet(1,1)&
&       *dgmetds(3,3)))
       cm(11,10,3,3)=5*gmet(2,3)**3*dgmetds(1,2)+7.5d0*gmet(1,3)*gmet(2,3)&
&       **2*dgmetds(2,2)-1.5d0*gmet(1,3)*gmet(2,2)*gmet(3,3)*dgmetds(2,2)&
&       -3*gmet(2,3)*gmet(3,3)*(gmet(2,2)*dgmetds(1,2)+gmet(1,2)*dgmetds(2,2))
       cm(12,10,3,3)=(gmet(2,3)*gmet(3,3)*(144*gmet(3,3)*dgmetds(1,2)&
&       +432*gmet(1,3)*dgmetds(2,3)-72*gmet(1,2)*dgmetds(3,3))+gmet(2,3)&
&       **2*(216*gmet(3,3)*dgmetds(1,3)+180*gmet(1,3)*dgmetds(3,3))+gmet(3,3)&
&       *(gmet(3,3)*(72*gmet(1,3)*dgmetds(2,2)-144*gmet(1,2)*dgmetds(2,3))&
&       +gmet(2,2)*(-72*gmet(3,3)*dgmetds(1,3)-36*gmet(1,3)*dgmetds(3,3))))&
&       /24.d0
       cm(13,10,3,3)=(180*gmet(2,3)**3*dgmetds(1,3)+gmet(2,3)*gmet(3,3)&
&       *(-108*gmet(2,2)*dgmetds(1,3)+324*gmet(1,3)*dgmetds(2,2)-216*gmet(1,2)&
&       *dgmetds(2,3))+gmet(2,3)**2*(324*gmet(3,3)*dgmetds(1,2)+540*gmet(1,3)&
&       *dgmetds(2,3))-108*gmet(3,3)*(gmet(1,2)*gmet(3,3)*dgmetds(2,2)&
&       +gmet(2,2)*(gmet(3,3)*dgmetds(1,2)+gmet(1,3)*dgmetds(2,3))))&
&       /36.d0
       cm(14,10,3,3)=gmet(3,3)*(2*gmet(3,3)**2*dgmetds(1,2)+9*gmet(1,3)&
&       *gmet(2,3)*dgmetds(3,3)+gmet(3,3)*(6*gmet(2,3)*dgmetds(1,3)+6*gmet(1,3)&
&       *dgmetds(2,3)-3*gmet(1,2)*dgmetds(3,3)))
       cm(15,10,3,3)=gmet(3,3)**2*(2*gmet(3,3)*dgmetds(1,3)+3*gmet(1,3)&
&       *dgmetds(3,3))
       cm(16,10,3,3)=gmet(2,3)*(2.5d0*gmet(2,3)**2-1.5d0*gmet(2,2)*gmet(3,3))&
&       *dgmetds(2,2)
       cm(17,10,3,3)=9*gmet(2,3)**2*gmet(3,3)*dgmetds(2,3)-3*gmet(2,2)&
&       *gmet(3,3)**2*dgmetds(2,3)+2.5d0*gmet(2,3)**3*dgmetds(3,3)+gmet(2,3)&
&       *gmet(3,3)*(3*gmet(3,3)*dgmetds(2,2)-1.5d0*gmet(2,2)*dgmetds(3,3))
       cm(18,10,3,3)=4.5d0*gmet(2,3)**2*gmet(3,3)*dgmetds(2,2)-1.5d0*gmet(2,2)&
&       *gmet(3,3)**2*dgmetds(2,2)+5*gmet(2,3)**3*dgmetds(2,3)-3*gmet(2,2)&
&       *gmet(2,3)*gmet(3,3)*dgmetds(2,3)
       cm(19,10,3,3)=gmet(3,3)*(1*gmet(3,3)**2*dgmetds(2,2)+4.5d0*gmet(2,3)&
&       **2*dgmetds(3,3)+gmet(3,3)*(6*gmet(2,3)*dgmetds(2,3)-1.5d0*gmet(2,2)&
&       *dgmetds(3,3)))
       cm(20,10,3,3)=gmet(3,3)**2*(2*gmet(3,3)*dgmetds(2,3)+3*gmet(2,3)&
&       *dgmetds(3,3))
       cm(21,10,3,3)=gmet(3,3)**3*dgmetds(3,3)


     end if
   end if
 end if !cm_set==0

!
!this is the part of the routine which gets executed on each call
!
 rankin=rank
 rankout=rank
 if(iterm==1) rankout=rankout+2
 if(iterm==3) rankin=rankin+2
 limitin=(rankin+1)*(rankin+2)/2
 limitout=(rankout+1)*(rankout+2)/2
!matrix-vector multiplication is small and probably best left written
!out be hand rather than using a LAPACK call
 do jj=1,limitout
   bb(:,jj)=0.d0
   do ii=1,limitin
     bb(:,jj)=bb(:,jj)+aa(:,ii)*cm(ii,jj,iterm,rank)
   end do
 end do

end subroutine metstr
!!***

end module m_metstr
!!***

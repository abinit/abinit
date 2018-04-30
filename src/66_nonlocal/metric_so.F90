!{\src2tex{textfont=tt}}
!!****f* ABINIT/metric_so
!! NAME
!! metric_so
!!
!! FUNCTION
!! Computes Pauli matrices and antisymmetric tensor needed for
!! spin-orbit.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GZ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!              (bohr**-1)
!!
!! OUTPUT
!!  amet(2,3,3,2,2)=the antisymmetric tensor A(Re/Im,y,y'',s,s'')
!!  pauli(2,2,2,3)=Pauli matrixes
!!
!! NOTES
!! (the preprocessing based on CPP does not allow isolated quotes,
!!  so that they have been doubled in the following latex equations)
!!{{\ \begin{eqnarray}
!! $2 Pauli(s,s'',1) & = & 0 & 1 \nonumber
!!                   &   & 1 & 0 \nonumber
!!  2 Pauli(s,s'',2) & = & 0 & -i \nonumber
!!                   &   & i &  0 \nonumber
!!  2 Pauli(s,s'',3) & = & 1 &  0 \nonumber
!!                   &   & 0 & -1
!! \end{eqnarray} }}
!!{{\ \begin{eqnarray}
!!    $Amet(y,y'',s,s'') & = & -i Pauli(s,s'',n) E(n,m,m'') Gprimd(m,y) Gprimd(m'',y'')
!! \end{eqnarray} }}
!!
!! E(n,m,m''):   full antisymetric tensor
!! s,s'':        spin indices (1..2)
!! y,y'',m,m'',n: metric indices (1..3)
!! a,b:         strain indices (1..3)
!! Amet and Pauli are complex
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


subroutine metric_so(amet,gprimd,pauli)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'metric_so'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(out) :: amet(2,3,3,2,2),pauli(2,2,2,3)

!Local variables-------------------------------
!scalars
 integer :: iy1,iy2,m1,m2,n
!arrays
 real(dp) :: buffer1(3,3,2,2) !,buffer2(3,3,3,3,2,2)

! **********************************************************************

!Fill in Pauli matrices and make them spin matrices:
!Pauli(Re/Im,up,down,coord):
 pauli(:,:,:,:)=0.d0
 pauli(1,1,2,1)= 1.d0;pauli(1,2,1,1)= 1.d0
 pauli(2,1,2,2)=-1.d0;pauli(2,2,1,2)= 1.d0
 pauli(1,1,1,3)= 1.d0;pauli(1,2,2,3)=-1.d0
 pauli(:,:,:,:)= 0.5d0*pauli(:,:,:,:)

!Construct the antisymmetric tensor:
 amet(:,:,:,:,:)=0.d0
 do iy2=1,3
   do iy1=1,3
     do n=1,3
       m1=mod(n ,3)+1    !  n,m1,m2 is an even permutation
       m2=mod(m1,3)+1
       amet(1:2,iy1,iy2,1:2,1:2) = amet(:,iy1,iy2,:,:) &
&       + pauli(:,:,:,n) &
&       *(gprimd(m1,iy1)*gprimd(m2,iy2) &
&       -gprimd(m2,iy1)*gprimd(m1,iy2))
     end do
   end do
 end do
!amet= -i amet
 buffer1(:,:,:,:)=-amet(1,:,:,:,:)
 amet(1,:,:,:,:) = amet(2,:,:,:,:)
 amet(2,:,:,:,:) =buffer1(:,:,:,:)

!DEBUG
!!Eventually construct the gradients of Amet wrt strains:
!! DAmet(y,y'',a,b,s,s'')= d[Amet(y,y'',s,s'')]/d[strain(a,b)]
!!                        -i Pauli(s,s'',n)*
!!                          ( E(n,a,m'')*Gprimd(b,y )*Gprimd(m'',y'')
!!                           +E(n,m,a )*Gprimd(b,y'')*Gprimd(m ,y ) )
!damet(:,:,:,:,:,:,:)=0.d0
!do ib=1,3
!do ia=1,3
!m1=mod(ia,3)+1    !  ia,m1,m2 is an even permutation
!m2=mod(m1,3)+1
!do iy2=1,3
!do iy1=1,3
!damet(:,iy1,iy2,ia,ib,:,:) = damet(:,iy1,iy2,ia,ib,:,:) &
!&        + (pauli(:,:,:,m2)*gprimd(m1,iy2) &
!-pauli(:,:,:,m1)*gprimd(m2,iy2))*gprimd(ib,iy1) &
!&        + (pauli(:,:,:,m1)*gprimd(m2,iy1) &
!-pauli(:,:,:,m2)*gprimd(m1,iy1))*gprimd(ib,iy2)
!end do
!end do
!end do
!end do
!! damet= i damet
!buffer2(:,:,:,:,:,:)= damet(1,:,:,:,:,:,:)
!damet(1,:,:,:,:,:,:)= -damet(2,:,:,:,:,:,:)
!damet(2,:,:,:,:,:,:)=buffer2(:,:,:,:,:,:)
!! Symetrize damet(:,:,:,a,b,:,:)
!damet(:,:,:,1,2,:,:)=0.5d0*(damet(:,:,:,1,2,:,:)+damet(:,:,:,2,1,:,:))
!damet(:,:,:,1,3,:,:)=0.5d0*(damet(:,:,:,1,3,:,:)+damet(:,:,:,3,1,:,:))
!damet(:,:,:,2,3,:,:)=0.5d0*(damet(:,:,:,2,3,:,:)+damet(:,:,:,3,2,:,:))
!damet(:,:,:,2,1,:,:)=damet(:,:,:,1,2,:,:)
!damet(:,:,:,3,1,:,:)=damet(:,:,:,1,3,:,:)
!damet(:,:,:,3,2,:,:)=damet(:,:,:,2,3,:,:)
!ENDDEBUG

end subroutine metric_so
!!***

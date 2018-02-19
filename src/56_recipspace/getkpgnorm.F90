!{\src2tex{textfont=tt}}
!!****f* ABINIT/getkpgnorm
!! NAME
!! getkpgnorm
!!
!! FUNCTION
!!  compute the norms of the k+G vectors
!!
!! COPYRIGHT
!! Copyright (C) 2003-2018 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gprimd(3,3)=metric tensor
!!  kg_k(3,npw_k)= G vectors, in reduced coordinates
!!  kpt(3)=k vector, in reduced coordinates
!!  npw_k=size of the G-vector set
!!
!! OUTPUT
!!  kpgnorm(npw_k)=norms of the k+G vectors
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_cut3d,partial_dos_fractions
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine getkpgnorm(gprimd,kpt,kg_k,kpgnorm,npw_k)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getkpgnorm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(in) :: gprimd(3,3),kpt(3)
 real(dp),intent(out) :: kpgnorm(npw_k)

!Local variables-------------------------------
!character(len=500) :: message
!scalars
 integer :: ipw
 real(dp) :: g11,g12,g13,g21,g22,g23,g31,g32,g33,k1,k2,k3,kpg1,kpg2,kpg3,rr,xx
 real(dp) :: yy,zz

! *************************************************************************

 k1=kpt(1) ; k2=kpt(2) ; k3=kpt(3)
 g11=gprimd(1,1)
 g12=gprimd(1,2)
 g13=gprimd(1,3)
 g21=gprimd(2,1)
 g22=gprimd(2,2)
 g23=gprimd(2,3)
 g31=gprimd(3,1)
 g32=gprimd(3,2)
 g33=gprimd(3,3)

!Loop over all k+G
 do ipw=1,npw_k

!  Load k+G
   kpg1=k1+dble(kg_k(1,ipw))
   kpg2=k2+dble(kg_k(2,ipw))
   kpg3=k3+dble(kg_k(3,ipw))

!  Calculate module of k+G
   xx=g11*kpg1+g12*kpg2+g13*kpg3
   yy=g21*kpg1+g22*kpg2+g23*kpg3
   zz=g31*kpg1+g32*kpg2+g33*kpg3
   rr=sqrt(xx**2+yy**2+zz**2)
   kpgnorm(ipw) = rr

 end do ! ipw

end subroutine getkpgnorm
!!***

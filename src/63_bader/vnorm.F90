!{\src2tex{textfont=tt}}
!!****f* ABINIT/vnorm
!! NAME
!! vnorm
!!
!! FUNCTION
!! Default declarations, and interfaces for the aim.f utility.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2017 ABINIT group (PCasek,FF,XG,TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! vector norm ->dir==1: vector in reduced coordinates
!!               dir==0: vector in cartes. coordinates
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! vv = on entry, vector to normalized
!!    = on exit, normalized vector
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

function vnorm(vv,dir)

 use defs_basis
 use defs_aimprom
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vnorm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dir
 real(dp) :: vnorm
!arrays
 real(dp),intent(in) :: vv(3)

!Local variables-------------------------------
!scalars
 integer :: ii
!arrays
 real(dp) :: vt(3)

! *************************************************************************

 vnorm=zero
 if (dir==1) then
   do ii=1,3
     vt(ii)=rprimd(ii,1)*vv(1)+rprimd(ii,2)*vv(2)+rprimd(ii,3)*vv(3)
     vnorm=vnorm+vt(ii)*vt(ii)
   end do
 elseif (dir==0) then
   do ii=1,3
     vnorm=vnorm+vv(ii)*vv(ii)
   end do
 else
   MSG_ERROR('vnorm calcul')
 end if
 vnorm=sqrt(vnorm)
end function vnorm
!!***


!!****f* ABINIT/vec_prod
!! NAME
!! vec_prod
!!
!! FUNCTION
!! Vector product
!!
!! COPYRIGHT
!! Copyright (C) 2007-2017 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  vv,uu = vectors to compute vector product
!! OUTPUT
!!  (return the value of the vector product)
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function vec_prod(uu,vv)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vec_prod'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp) :: vec_prod(3)
 real(dp),intent(in) :: uu(3),vv(3)

!Local variables-------------------------------

! *************************************************************************

 vec_prod(1)=uu(2)*vv(3)-vv(2)*uu(3)
 vec_prod(2)=uu(3)*vv(1)-vv(3)*uu(1)
 vec_prod(3)=uu(1)*vv(2)-vv(1)*uu(2)

end function vec_prod
!!***


!!****f* ABINIT/mprod
!! NAME
!! mprod
!!
!! FUNCTION
!! Matrix multiplication cc=aa*bb
!!
!! PARENTS
!!      vnorm
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mprod(aa,bb,cc)

 use m_profiling_abi

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mprod'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: aa(3,3),bb(3,3)
 real(dp),intent(out) :: cc(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,kk

! *************************************************************************

 do ii=1,3
   do jj=1,3
     cc(ii,jj)=0._dp
     do kk=1,3
       cc(ii,jj)=cc(ii,jj)+aa(ii,kk)*bb(kk,jj)
     end do
   end do
 end do
end subroutine mprod
!!***


!!****f* ABINIT/bschg1
!! NAME
!! bschg1
!!
!! FUNCTION
!! bschg1: Vector transformation of coordinates
!!
!! PARENTS
!!      critics,initaim,integrho,vgh_rho
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine bschg1(vv,dir)

 use defs_basis
 use defs_aimprom
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bschg1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dir
!arrays
 real(dp),intent(inout) :: vv(3)

!Local variables ------------------------------
!scalars
 integer :: ii
!arrays
 real(dp) :: vt(3)

! *********************************************************************

 if (dir==1) then
   do ii=1,3
     vt(ii)=rprimd(ii,1)*vv(1)+rprimd(ii,2)*vv(2)+rprimd(ii,3)*vv(3)
   end do
 elseif (dir==-1) then
   do ii=1,3
     vt(ii)=ivrprim(ii,1)*vv(1)+ivrprim(ii,2)*vv(2)+ivrprim(ii,3)*vv(3)
   end do
 elseif (dir==2) then
   do ii=1,3
     vt(ii)=trivrp(ii,1)*vv(1)+trivrp(ii,2)*vv(2)+trivrp(ii,3)*vv(3)
   end do
 else
   MSG_ERROR('Transformation of coordinates')
 end if
 vv(:)=vt(:)
end subroutine bschg1
!!***

!!****f* ABINIT/bschg2
!! NAME
!! bschg2
!!
!! FUNCTION
!! bschg2: Matrix transformation of coordinates
!!
!! PARENTS
!!      vgh_rho
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine bschg2(aa,dir)

 use defs_basis
 use defs_aimprom
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bschg2'
 use interfaces_63_bader, except_this_one => bschg2
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dir
!arrays
 real(dp),intent(inout) :: aa(3,3)

!Local variables ------------------------------
!arrays
 real(dp) :: bb(3,3)

! *********************************************************************

 if (dir==1) then
   call mprod(aa,ivrprim,bb)
   call mprod(rprimd,bb,aa)
 elseif (dir==2) then
   call mprod(aa,ivrprim,bb)
   call mprod(trivrp,bb,aa)
 elseif (dir==-1) then
   call mprod(aa,rprimd,bb)
   call mprod(ivrprim,bb,aa)
 else
   MSG_ERROR("transformation of coordinates")
 end if
end subroutine bschg2
!!***

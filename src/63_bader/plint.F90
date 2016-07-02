!{\src2tex{textfont=tt}}
!!****f* ABINIT/plint
!! NAME
!! plint
!!
!! FUNCTION
!! This simple routine gives the profile of the density
!! integrated in xy plane belong the z-axes (it works only
!! for orthogonal coordinates at present - it is better to use cut3d)
!! integration in plane - with equilateral triangles (not really
!! finished and not tested!)
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (this routine works on the data in the aimprom module)
!!
!! OUTPUT
!!  (this routine works on the data in the aimprom module)
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      drvaim
!!
!! CHILDREN
!!      vgh_rho
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine plint()

 use m_profiling_abi

 use defs_basis
 use defs_aimprom

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'plint'
 use interfaces_63_bader, except_this_one => plint
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!Local variables ------------------------------
!scalars
 integer,parameter :: nd=150,ng=300
 integer :: cod,iat,ii,ipos,jj,kk,nn
 real(dp) :: dd,ee,ff,gg,hh,igr,rho,ss
 logical :: prep
!arrays
 real(dp) :: grho(3),hrho(3,3),vv(3),xl(nd+1),xs(nd)
 real(dp),allocatable :: uu(:)

! *********************************************************************

 ff=rprimd(1,1)/nd
 ss=2._dp/sqrt(3._dp)*rprimd(2,2)/rprimd(1,1)*nd
 nn=int(ss)
 gg=sqrt(3._dp)/2.*ff
 hh=rprimd(2,2)-nn/nd*sqrt(3._dp)/2.*rprimd(1,1)
 ee=hh/sqrt(3._dp)
 hh=hh/2.
 ss=sqrt(3._dp)*ff*ff/4.
 dd=ee*ff/2.

 do ii=1,nd
   xl(ii)=ii*ff
   xs(ii)=ff/2.+ii*ff
 end do
 xl(nd+1)=rprimd(1,1)

 ABI_ALLOCATE(uu,(nn+3))

 uu(1)=0._dp
 uu(nn+3)=rprimd(2,2)
 do ii=2,nn+2
   uu(ii)=hh+(ii-1)*gg
 end do
 igr=0._dp
 prep=.true.
 do kk=1,ng
   igr=0._dp
   vv(3)=(kk-1)*rprimd(3,3)/ng
   do ii=1,nn+3
     vv(2)=uu(ii)
     do jj=1,nd
       if (prep) then
         vv(1)=xl(jj)
         prep=.false.
       else
         vv(1)=xs(jj)
         prep=.true.
       end if
       call vgh_rho(vv,rho,grho,hrho,dd,iat,ipos,cod)
       if ((ii==1).or.(ii==nn+3)) then
         igr=igr+dd*rho
       elseif ((ii==2).or.(ii==nn+2)) then
         igr=igr+(dd+ss)*rho
       else
         igr=igr+ss*2*rho
       end if
     end do
   end do
   write(untp,'(2E16.8)') vv(3), igr
 end do
 ABI_DEALLOCATE(uu)
end subroutine plint
!!***

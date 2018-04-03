!{\src2tex{textfont=tt}}
!!****f* ABINIT/vhtnzc
!! NAME
!! vhtnzc
!!
!! FUNCTION
!! Compute VHartree(tild[n_Z+n_core]) from input ncore
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (GJ, FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mesh=dimension of radial mesh
!!  nc= core density (to be pseudized)
!!  rad(mesh)=radial mesh
!!  rc=cut-off radius
!!  znucl=nuclear number of atom as specified in psp file
!!
!! OUTPUT
!!  vhtnzc(mesh) = hartree potential induced by density tild[n_Z+n_core] (pseudo core density + nucleus)
!!
!! PARENTS
!!      psp6cc
!!
!! CHILDREN
!!      ctrap
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine vhtnzc(nc,rc,vh_tnzc,mesh,rad,znucl)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vhtnzc'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mesh
 real(dp),intent(in) :: znucl
 real(dp),intent(in) :: rc
!arrays
 real(dp),intent(in) :: nc(mesh),rad(mesh)
 real(dp),intent(out) :: vh_tnzc(mesh)

!Local variables-------------------------------
!scalars
 integer :: ir,nc1
 real(dp) :: gnorm,rc1,step,yp1,yp2,yp3
!arrays
 real(dp),allocatable :: den1(:),den2(:),den3(:),den4(:),nzc(:),rvhn(:),shapefunc(:)

! *************************************************************************

 rc1=rc/four

 step=log(rad(2)/rad(1))
 nc1=int(log(rc1/rad(1))/step)+1
 rc1=rad(nc1)

 ABI_ALLOCATE(shapefunc,(mesh))
 shapefunc(1)=one
 shapefunc(2:nc1)=(sin(pi*rad(2:nc1)/rc1)/(pi*rad(2:nc1)/rc1))**2
 if (nc1<mesh) shapefunc(nc1+1:mesh)=zero

 ABI_ALLOCATE(den1,(mesh))
 den1(1:mesh)=four_pi*shapefunc(1:mesh)*rad(1:mesh)**3
 call ctrap(mesh,den1,step,gnorm)
 gnorm =one/gnorm
 ABI_DEALLOCATE(den1)

 ABI_ALLOCATE(nzc,(mesh))
 nzc(1:mesh)=four*pi*nc(1:mesh)*rad(1:mesh)**2-four_pi*shapefunc(1:mesh)*rad(1:mesh)**2*znucl*gnorm
 ABI_DEALLOCATE(shapefunc)

 ABI_ALLOCATE(rvhn,(mesh))
 rvhn(1)=zero

 ABI_ALLOCATE(den1,(mesh))
 ABI_ALLOCATE(den2,(mesh))
 ABI_ALLOCATE(den3,(mesh))
 ABI_ALLOCATE(den4,(mesh))

 den1(1)=zero;den2(1)=zero
 do ir=2,mesh
   den1(ir)= rad(ir)*nzc(ir)
   den2(ir)= den1(ir)/rad(ir)
 end do

!For first few points do stupid integral
 den3(1)=zero;den4(1)=zero
 do ir=2,mesh
   call ctrap(ir,den1(1:ir),step,den3(ir))
   call ctrap(ir,den2(1:ir),step,den4(ir))
 end do

 do ir=1,mesh
   rvhn(ir)=den3(ir)+rad(ir)*(den4(mesh)-den4(ir))
 end do

 ABI_DEALLOCATE(den1)
 ABI_DEALLOCATE(den2)
 ABI_DEALLOCATE(den3)
 ABI_DEALLOCATE(den4)

 vh_tnzc(2:mesh)=rvhn(2:mesh)/rad(2:mesh)
 yp2=(vh_tnzc(3)-vh_tnzc(2))/(rad(3)-rad(2))
 yp3=(vh_tnzc(4)-vh_tnzc(3))/(rad(4)-rad(3))
 yp1=yp2+(yp2-yp3)*rad(2)/(rad(3)-rad(2))
 vh_tnzc(1)=vh_tnzc(2)-(yp1+yp2)*rad(2)

 ABI_DEALLOCATE(nzc)
 ABI_DEALLOCATE(rvhn)

end subroutine vhtnzc
!!***

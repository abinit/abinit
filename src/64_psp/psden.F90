!{\src2tex{textfont=tt}}
!!****f* ABINIT/psden
!! NAME
!! psden
!!
!! FUNCTION
!! Calculate a pseudo-density from an original density on a radial grid (regular or logarithmic)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (GJ,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ilog=1 if grid is logarithmic, else 0
!!  mesh= dimension of nc
!!  nc(mesh)= density to be pseudized
!!  rc= cut-off radius
!!  rad(mesh) = radial mesh
!!
!! OUTPUT
!!  ff(mesh)= pseudized density
!!
!!SIDE EFFECTS
!!  Optional:
!!    ff1(mesh)= 1st derivative of pseudo density (only r<rc modified)
!!    ff2(mesh)= 2nd derivative of pseudo density (only r<rc modified)
!!
!! NOTES
!!    ff=exp(-(a+b.r^2+c.r^4))
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

subroutine psden(ilog,ff,mesh,nc,rc,rad,ff1,ff2)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_numeric_tools, only : ctrap

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psden'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ilog,mesh
 real(dp),intent(in) :: rc
!arrays
 real(dp),intent(in) :: nc(mesh),rad(mesh)
 real(dp),intent(out) :: ff(mesh)
 real(dp),intent(inout),optional :: ff1(mesh),ff2(mesh)

!Local variables-------------------------------
!scalars
 integer :: ii,nc1
 real(dp) :: aa,aa1,aa2,bb,cc,c1,c3,f0,f0p,norm1,norm2,rc1,step
!arrays
 real(dp),allocatable :: fpir(:),gg(:)

! *************************************************************************

 rc1=rc/four

 ABI_ALLOCATE(fpir,(mesh))
 fpir(1:mesh)=four_pi*rad(1:mesh)**2
 if (ilog==1) fpir(1:mesh)=fpir(1:mesh)*rad(1:mesh)

 if (ilog==0) then
   step=rad(2)-rad(1)
   nc1=int(rc1/step)+1
   rc1=(nc1-1)*step
 else if (ilog==1) then
   step=log(rad(2)/rad(1))
   nc1=int(log(rc1/rad(1))/step)+1
   rc1=rad(nc1)
 end if
 ff(1:nc1)=nc(1:nc1)*fpir(1:nc1)
 call ctrap(nc1,ff(1:nc1),step,c3)
 if (ilog==1) c3=c3+half*ff(1)
 f0=nc(nc1);c1=-log(f0)
 f0p=half*(nc(nc1+1)-nc(nc1-1))/step

 ii=0;aa1=zero;norm1=c3+one
 do while (norm1>c3.and.ii<100)
   ii=ii+1;aa1=aa1+one
   aa=c1-aa1*rc1**4+rc1*(f0p/f0+four*aa1*rc1**3)*half
   bb=-half*(f0p/f0+four*aa1*rc1**3)/rc1
   ff(1:nc1)=fpir(1:nc1)*exp(-aa-bb*rad(1:nc1)**2-aa1*rad(1:nc1)**4)
   call ctrap(nc1,ff(1:nc1),step,norm1)
   if (ilog==1) norm1=norm1+half*ff(1)
 end do
 if (ii==100) then
   MSG_ERROR('Big pb 1 in psden !')
 end if

 ii=0;aa2=zero;norm2=c3-one
 do while (norm2<c3.and.ii<100)
   ii=ii+1;aa2=aa2-one
   aa=c1-aa2*rc1**4+rc1*(f0p/f0+four*aa2*rc1**3)*half
   bb=-half*(f0p/f0+four*aa2*rc1**3)/rc1
   ff(1:nc1)=fpir(1:nc1)*exp(-aa-bb*rad(1:nc1)**2-aa2*rad(1:nc1)**4)
   call ctrap(nc1,ff(1:nc1),step,norm2)
   if (ilog==1) norm2=norm2+half*ff(1)
 end do
 if (ii==100) then
   MSG_ERROR('Big pb 2 in psden !')
 end if

 do while (abs(norm2-c3)>tol10)

   cc=(aa1+aa2)*half
   aa=c1-cc*rc1**4+rc1*(f0p/f0+four*cc*rc1**3)*half
   bb=-half*(f0p/f0+four*cc*rc1**3)/rc1
   ff(1:nc1)=fpir(1:nc1)*exp(-aa-bb*rad(1:nc1)**2-cc*rad(1:nc1)**4)
   call ctrap (nc1,ff(1:nc1),step,norm2)
   if (ilog==1) norm2=norm2+half*ff(1)
   if ((norm1-c3)*(norm2-c3)>zero) then
     aa1=cc
     norm1=norm2
   else
     aa2=cc
   end if

 end do ! while

 ff(1)=exp(-aa);if (ilog==1) ff(1)=ff(1)*exp(-bb*rad(1)**2-cc*rad(1)**4)
 ff(2:nc1)=ff(2:nc1)/fpir(2:nc1)
 if (nc1<mesh) ff(nc1+1:mesh)=nc(nc1+1:mesh)
 if (present(ff1)) ff1(1:nc1)=-(two*bb*rad(1:nc1)+four*cc*rad(1:nc1)**3)*ff(1:nc1)
 if (present(ff2)) ff2(1:nc1)=-(two*bb+12.0_dp*cc*rad(1:nc1)**2)*ff(1:nc1) &
& +(two*bb*rad(1:nc1)+four*cc*rad(1:nc1)**3)**2*ff(1:nc1)

 ABI_ALLOCATE(gg,(mesh))
 gg(1:mesh)=fpir(1:mesh)*ff(1:mesh)
 call ctrap(mesh,gg(1:mesh),step,norm1)
 if (ilog==1) norm1=norm1+half*gg(1)
 write(std_out,*) 'psden: tild_nc integral= ',norm1
 ABI_DEALLOCATE(gg)

 ABI_DEALLOCATE(fpir)

end subroutine psden
!!***

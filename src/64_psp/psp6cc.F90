!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp6cc
!! NAME
!! psp6cc
!!
!! FUNCTION
!! Compute the core charge density, for use in the XC core
!! correction, following the function definition valid
!! for the format 6 of pseudopotentials.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (AF, GJ,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mmax=maximum number of points in real space grid in the psp file
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  rchrg=cut-off radius for the core density
!!  znucl=nuclear number of atom as specified in psp file
!!
!! OUTPUT
!!  xccc1d(n1xccc,6)= 1D core charge function and its five first derivatives
!!  Optional output:
!!    vh_tnzc(mmax) = Hartree potential induced by density tild_[n_Z+n_core]
!!                    (pseudized [n_Z+n_core], where n_Z=ions, n_core=core electrons)
!!                    using a simple pseudization scheme
!!
!! PARENTS
!!      psp6in
!!
!! CHILDREN
!!      psden,smooth,spline,splint,vhtnzc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp6cc(mmax,n1xccc,rchrg,xccc1d,znucl,&
&                 vh_tnzc) ! optional argument

 use defs_basis
 use m_splines
 use m_profiling_abi

 use m_numeric_tools,  only : smooth

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp6cc'
 use interfaces_64_psp, except_this_one => psp6cc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mmax,n1xccc
 real(dp),intent(in) :: rchrg,znucl
!arrays
 real(dp),intent(inout) :: xccc1d(n1xccc,6) !vz_i
 real(dp),intent(out),optional :: vh_tnzc(mmax)

!Local variables-------------------------------
!scalars
 integer :: i1xccc,irad
 real(dp) :: der1,dern
 character(len=500) :: errmsg
!arrays
 real(dp),allocatable :: ff(:),ff1(:),ff2(:),ff3(:),gg(:),gg1(:),gg2(:),gg3(:)
 real(dp),allocatable :: gg4(:),nc(:),rad(:),work(:),xx(:)

!**********************************************************************

 ABI_ALLOCATE(ff,(mmax))
 ABI_ALLOCATE(ff1,(mmax))
 ABI_ALLOCATE(ff2,(mmax))
 ABI_ALLOCATE(ff3,(mmax))
 ABI_ALLOCATE(rad,(mmax))
 ABI_ALLOCATE(gg,(n1xccc))
 ABI_ALLOCATE(gg1,(n1xccc))
 ABI_ALLOCATE(gg2,(n1xccc))
 ABI_ALLOCATE(gg3,(n1xccc))
 ABI_ALLOCATE(gg4,(n1xccc))
 ABI_ALLOCATE(work,(n1xccc))
 ABI_ALLOCATE(xx,(n1xccc))

!
!read from pp file the model core charge (ff) and first (ff1) and
!second (ff2) derivative on logarithmic mesh mmax; rad is the radial grid
!the input functions contain the 4pi factor, it must be rescaled.

 do irad=1,mmax
   read(tmp_unit,*, err=10, iomsg=errmsg) rad(irad),ff(irad),ff1(irad),ff2(irad)
   ff(irad)=ff(irad)/four_pi
   ff1(irad)=ff1(irad)/four_pi
   ff2(irad)=ff2(irad)/four_pi
 end do

!Optional output: VHartree(tild_[n_Z+n_core])
 if (present(vh_tnzc)) then
   ABI_ALLOCATE(nc,(mmax))
   nc=ff ! n_core
   call psden(1,ff,mmax,nc,rchrg,rad,ff1=ff1,ff2=ff2)
   call vhtnzc(ff,rchrg,vh_tnzc,mmax,rad,znucl)
   ABI_DEALLOCATE(nc)
 end if

 rad(1)=zero

!calculate third derivative ff3 on logarithmic grid
 der1=ff2(1)
 dern=ff2(mmax)
 call spline(rad,ff1,mmax,der1,dern,ff3)

!generate uniform mesh xx in the box cut by rchrg:

 do i1xccc=1,n1xccc
   xx(i1xccc)=(i1xccc-1)* rchrg/dble(n1xccc-1)
 end do

!
!now interpolate core charge and derivatives on the uniform grid
!
!core charge, input=ff,  output=gg
 call splint(mmax,rad,ff,ff2,n1xccc,xx,gg)

!first derivative input=ff1, output=gg1
 call splint(mmax,rad,ff1,ff3,n1xccc,xx,gg1)

!normalize gg1
 gg1(:)=gg1(:)*rchrg

!now calculate second to fourth derivative by forward differences
!to avoid numerical noise uses a smoothing function

 call smooth(gg1,n1xccc,10)

 gg2(n1xccc)=zero
 do i1xccc=1,n1xccc-1
   gg2(i1xccc)=(gg1(i1xccc+1)-gg1(i1xccc))*dble(n1xccc-1)
 end do

 call smooth(gg2,n1xccc,10)

 gg3(n1xccc)=zero
 do i1xccc=1,n1xccc-1
   gg3(i1xccc)=(gg2(i1xccc+1)-gg2(i1xccc))*dble(n1xccc-1)
 end do

 call smooth(gg3,n1xccc,10)

 gg4(n1xccc)=zero
 do i1xccc=1,n1xccc-1
   gg4(i1xccc)=(gg3(i1xccc+1)-gg3(i1xccc))*dble(n1xccc-1)
 end do

 call smooth(gg4,n1xccc,10)

!write on xcc1d
 xccc1d(:,1)=gg(:)
 xccc1d(:,2)=gg1(:)
 xccc1d(:,3)=gg2(:)
 xccc1d(:,4)=gg3(:)
 xccc1d(:,5)=gg4(:)

!WARNING : fifth derivative not yet computed
 xccc1d(:,6)=zero

!DEBUG
!note: the normalization condition is the following:
!4pi rchrg /dble(n1xccc-1) sum xx^2 xccc1d(:,1) = qchrg
!
!norm=zero
!do i1xccc=1,n1xccc
!norm = norm + four_pi*rchrg/dble(n1xccc-1)*&
!&             xx(i1xccc)**2*xccc1d(i1xccc,1)
!end do
!write(std_out,*) ' norm=',norm
!
!write(std_out,*)' psp1cc : output of core charge density and derivatives '
!write(std_out,*)'   xx          gg           gg1  '
!do i1xccc=1,n1xccc
!write(10, '(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,1),xccc1d(i1xccc,2)
!end do
!write(std_out,*)'   xx          gg2          gg3  '
!do i1xccc=1,n1xccc
!write(11, '(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,3),xccc1d(i1xccc,4)
!end do
!write(std_out,*)'   xx          gg4          gg5  '
!do i1xccc=1,n1xccc
!write(12, '(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,5),xccc1d(i1xccc,6)
!end do
!write(std_out,*)' psp1cc : debug done, stop '
!stop
!ENDDEBUG

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(ff1)
 ABI_DEALLOCATE(ff2)
 ABI_DEALLOCATE(ff3)
 ABI_DEALLOCATE(gg)
 ABI_DEALLOCATE(gg1)
 ABI_DEALLOCATE(gg2)
 ABI_DEALLOCATE(gg3)
 ABI_DEALLOCATE(gg4)
 ABI_DEALLOCATE(rad)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(xx)

 return

 ! Handle IO error
 10 continue
 MSG_ERROR(errmsg)

end subroutine psp6cc
!!***

!!****f* ABINIT/psden
!! NAME
!! psden
!!
!! FUNCTION
!! Calculate a pseudo-density from an original density on a radial grid (regular or logarithmic)
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

subroutine psden(ilog,ff,mesh,nc,rc,rad,ff1,ff2)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_numeric_tools, only : ctrap

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psden'
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

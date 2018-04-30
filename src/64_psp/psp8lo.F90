!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp8lo
!! NAME
!! psp8lo
!!
!! FUNCTION
!! Compute sine transform to transform from V(r) to q^2 V(q).
!! Computes integrals on linear grid interpolated from the linear input
!! grid with a spacing adjusted to ensure convergence at the maximum
!! wavevector using corrected trapezoidal integration.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DRH, DCA, XG, FrD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  amesh=spacing for linear radial atomic grid.
!!  mmax=number of radial r grid points
!!  mqgrid=number of grid points in q from 0 to qmax.
!!  qgrid(mqgrid)=q grid values (bohr**-1).
!!  rad(mmax)=r grid values (bohr).
!!  vloc(mmax)=V(r) on radial grid.
!!  zion=nominal valence charge of atom.
!!
!! OUTPUT
!!  epsatm=$ 4\pi\int[r^2 (V(r)+\frac{Zv}{r}dr]$.
!!{{\\ \begin{equation}
!!  q2vq(mqgrid)
!!   =q^2 V(q)
!!   = -\frac{Zv}{\pi}
!!     + q^2 4\pi\int[(\frac{\sin(2\pi q r)}{2\pi q r})(r^2 V(r)+r Zv)dr].
!!\end{equation} }}
!!  yp1,ypn=derivative of q^2 V(q) wrt q at q=0 and q=qmax
!!   (needed for spline fitter).
!!
!! PARENTS
!!      psp8in,psp9in
!!
!! CHILDREN
!!      ctrap,splfit,spline
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp8lo(amesh,epsatm,mmax,mqgrid,qgrid,q2vq,rad,vloc,yp1,ypn,zion)

 use defs_basis
 use m_splines
 use m_profiling_abi

 use m_numeric_tools, only : ctrap

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp8lo'
!End of the abilint section

 implicit none

!Arguments----------------------------------------------------------
!scalars
 integer,intent(in) :: mmax,mqgrid
 real(dp),intent(in) :: amesh,zion
 real(dp),intent(out) :: epsatm,yp1,ypn
!arrays
 real(dp),intent(in) :: qgrid(mqgrid),rad(mmax),vloc(mmax)
 real(dp),intent(out) :: q2vq(mqgrid)

!Local variables-------------------------------
!Following parameter controls accuracy of Fourier transform based on qmax
!and represents the minimun number of integration points in one period.
!scalars
 integer,parameter :: NPT_IN_2PI=200
 integer :: ider,iq,ir,irmu,irn,mesh_mult,mmax_new
 real(dp) :: amesh_new,arg,fp1,fpn,qmesh,result,ztor1
!arrays
 real(dp),allocatable :: rad_new(:),rvlpz(:),rvlpz_new(:),sprvlpz(:,:),work(:)

! *************************************************************************

 ABI_ALLOCATE(work,(mmax))
 ABI_ALLOCATE(rvlpz,(mmax))

!Do q=0 separately (compute epsatm)
 ztor1=(zion/2.0d0+rad(1)*vloc(1)/3.d0)*rad(1)**2
!Set up integrand for q=0: $ \int[r^2 (V(r)+\frac{Zv}{r}) dr]$
 do ir=1,mmax
   rvlpz(ir)=rad(ir)*vloc(ir)+zion
   work(ir)=rad(ir)*rvlpz(ir)
 end do

!Do integral from zero to r(max)
 call ctrap(mmax,work,amesh,result)

 epsatm=4.d0*pi*result
 q2vq(1)=-zion/pi

!Find r mesh spacing necessary for accurate integration at qmax
 amesh_new=2.d0*pi/(NPT_IN_2PI*qgrid(mqgrid))

!Choose submultiple of input mesh
 mesh_mult=int(amesh/amesh_new) + 1
 mmax_new=mesh_mult*(mmax-1)+1
 amesh_new=amesh/dble(mesh_mult)

 ABI_ALLOCATE(rad_new,(mmax_new))
 ABI_ALLOCATE(rvlpz_new,(mmax_new))

 if(mesh_mult==1) then
   rad_new(:)=rad(:)
   rvlpz_new(:)=rvlpz(:)
 else
!  Set up spline and interpolate to finer mesh.
!  First, compute derivatives at end points
   fp1=(-50.d0*rvlpz(1)+96.d0*rvlpz(2)-72.d0*rvlpz(3)+32.d0*rvlpz(4)&
&   -6.d0*rvlpz(5))/(24.d0*amesh)
   fpn=(6.d0*rvlpz(mmax-4)-32.d0*rvlpz(mmax-3)+72.d0*rvlpz(mmax-2)&
&   -96.d0*rvlpz(mmax-1)+50.d0*rvlpz(mmax))/(24.d0*amesh)
   ABI_ALLOCATE(sprvlpz,(mmax,2))
   work(:)=zero

!  Spline fit
   call spline(rad, rvlpz,mmax,fp1,fpn,sprvlpz(:,2))
   sprvlpz(:,1)=rvlpz(:)

!  Set up new radial mesh
   irn=1
   do ir=1,mmax-1
     do irmu=0,mesh_mult-1
       rad_new(irn)=rad(ir)+dble(irmu)*amesh_new
       irn=irn+1
     end do
   end do
   rad_new(mmax_new)=rad(mmax)

   ider=0
   call splfit(rad,work,sprvlpz,ider,rad_new,rvlpz_new,mmax,mmax_new)

   ABI_DEALLOCATE(sprvlpz)
   ABI_DEALLOCATE(work)
   ABI_ALLOCATE(work,(mmax_new))
 end if

!Loop over q values
 do iq=2,mqgrid
   arg=2.d0*pi*qgrid(iq)

!  Set up integrand
   do  ir=1,mmax_new
     work(ir)=sin(arg*rad_new(ir))*rvlpz_new(ir)
   end do

!  Do integral from zero to rad(mmax)
   call ctrap(mmax_new,work,amesh_new,result)

!  Store q^2 v(q)
   q2vq(iq)=q2vq(1)+2.d0*qgrid(iq)*result

 end do

!Compute derivatives of q^2 v(q) at ends of interval
 qmesh=qgrid(2)-qgrid(1)
 yp1=(-50.d0*q2vq(1)+96.d0*q2vq(2)-72.d0*q2vq(3)+32.d0*q2vq(4)&
& -6.d0*q2vq(5))/(24.d0*qmesh)
 ypn=(6.d0*q2vq(mqgrid-4)-32.d0*q2vq(mqgrid-3)+72.d0*q2vq(mqgrid-2)&
& -96.d0*q2vq(mqgrid-1)+50.d0*q2vq(mqgrid))/(24.d0*qmesh)

 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(rad_new)
 ABI_DEALLOCATE(rvlpz_new)
 ABI_DEALLOCATE(rvlpz)

end subroutine psp8lo
!!***

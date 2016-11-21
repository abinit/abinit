!{\src2tex{textfont=tt}}
!!****f* ABINIT/mlwfovlp_radial
!! NAME
!! mlwfovlp_radial
!!
!! FUNCTION
!! Calculates the radial part of the initial functions given to Wannier90 
!! as an starting point for the minimization.
!! The trial functions are a set of solutions to the radial part of the hydrogenic
!! Schrodinger equation as it is explained in Table 3.3 of the Wannier90 user guide.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2016 ABINIT group (trangel,drh)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  alpha= Z/a = zona
!!  lmax= maximum value of l
!!  rvalue= integer defining the choice for radial functions R(r).
!!   It can take values from 1-3. 
!!   It is associted to the radial part of the hydrogenic Schrodinger equation for l=0,
!!   See the manual of Wannier90 for more information. (www.wannier.org)
!!  xx= scalar number used to calculate the spherical bessel function. J_il(xx)
!!
!! OUTPUT
!!  mlwfovlp_radial= radial part for initial projections used to construct MLWF
!!
!! SIDE EFFECTS
!!  None
!!
!! NOTES
!!  Calculates the radial part of the initial functions given as an initial
!!  guess by the user to construct the MLWF.
!!  
!! PARENTS
!!      mlwfovlp_proj
!!
!! CHILDREN
!!      besjm,simpson_int
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mlwfovlp_radial(alpha,lmax,lmax2,radial,rvalue,xx)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_numeric_tools,   only : simpson_int
 use m_special_funcs,   only : besjm

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mlwfovlp_radial'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmax,lmax2,rvalue
 real(dp),intent(in) :: alpha,xx
!arrays
 real(dp),intent(out) :: radial(lmax2)

!Local variables
!scalars
 integer :: ir,ll,lm,mesh,mm
 real(dp),parameter :: dx=0.015d0,rmax=10.d0,xmin=0.d0
 real(dp) :: aa,ftmp,gauss,rtmp,x
 character(len=500) :: message
!arrays
 real(dp),save :: dblefact(4)=(/1_dp,3_dp,15_dp,105_dp/)
 real(dp),allocatable :: aux(:),bes(:),cosr(:),func_r(:),r(:),rad_int(:)
 real(dp),allocatable :: sinr(:)

! *************************************************************************
 
!Radial functions in the form of hydrogenic orbitals as defined in the 
!wannier90 manual.
 if(( rvalue > 0 ).and.(rvalue < 4)) then

!  mesh
   mesh= nint((rmax - xmin ) / dx + 1)
   ABI_ALLOCATE( bes,(mesh))
   ABI_ALLOCATE(func_r,(mesh))
   ABI_ALLOCATE(r,(mesh))
   ABI_ALLOCATE(rad_int,(mesh))
   ABI_ALLOCATE( aux,(mesh))
   ABI_ALLOCATE(cosr,(mesh))
   ABI_ALLOCATE(sinr,(mesh))
   do ir=1, mesh
     x=xmin+DBLE(ir-1)*dx
     r(ir)=x
   end do   !ir

!  radial functions shown in table 3.3 of wannier90 manual
   if (rvalue==1) func_r(:) = 2.d0 * alpha**(3.d0/2.d0) * exp(-alpha*r(:))
   if (rvalue==2) func_r(:) = 1.d0/(2.d0*sqrt(2.d0))*alpha**(3.d0/2.d0) *&
&   (2.d0 - alpha*r(:))*exp(-alpha*r(:)/2.d0)
   if (rvalue==3) func_r(:) = sqrt(4.d0/27.d0)*alpha**(3.d0/2.d0)&
&   * (1.d0 - 2.d0*alpha*r(:)/3.d0 + 2.d0*alpha**2*r(:)**2/27.d0)&
&   * exp(-alpha * r(:)/3.d0)

!  compute spherical bessel functions
   cosr(:)=cos(xx*r(:))
   sinr(:)=sin(xx*r(:))
   lm=0
   do ll=0,lmax
     call besjm(xx,bes,cosr,ll,mesh,sinr,r)
     aux(:)=bes(:)*func_r(:)*r(:)
!    do ir=1,mesh
!    write(310,*) r(ir),bes(ir)
!    end do
     call simpson_int(mesh,dx,aux,rad_int)
     rtmp=rad_int(mesh)/mesh
     do mm=-ll,ll
       lm=lm+1
       radial(lm)=rtmp
     end do !mm
   end do !ll
   ABI_DEALLOCATE(bes)
   ABI_DEALLOCATE(func_r)
   ABI_DEALLOCATE(r)
   ABI_DEALLOCATE(aux)
   ABI_DEALLOCATE(rad_int)
   ABI_DEALLOCATE(cosr)
   ABI_DEALLOCATE(sinr)

!  Radial part in the form of Gaussian functions of a given width 
!  Taken by code of made by drh.
 elseif ( rvalue == 4) then
   aa=1._dp/alpha
   gauss=exp(-0.25_dp*(aa*xx)**2)
   lm=0
   do ll=0,lmax
     ftmp=(0.5_dp*pi)**(0.25_dp)*aa*sqrt(aa/dblefact(ll+1))*(aa*xx)**ll*gauss
     do mm=-ll,ll
       lm=lm+1
       radial(lm)=ftmp
     end do
   end do
 else ! rvalue < 0 of rvalue > 4
   write(message,'(a,i6,5a)')&
&   '  Radial function r=',rvalue,ch10,&
&   '  is not defined',ch10,&
&   '  Modify .win file',ch10
   MSG_BUG(message)
 end if !rvalue

end subroutine mlwfovlp_radial
!!***

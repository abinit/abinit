!{\src2tex{textfont=tt}}
!!****f* ABINIT/etheta
!! NAME
!! etheta
!!
!! FUNCTION
!! Computes the energy per unit cell and its first derivative
!! for a given angle theta. More precisely, computes only the part of
!! the energy that changes with theta.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2018 ABINIT  group (MVeithen,ISouza,JIniguez)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! bcut(ifor,idir) = branch cut of the ellipse associated with (ifor,idir)
!! chc = <C|H_0|C> where |C> is the wavefunction of the current band
!! detovc = determinant of the overlap matrix S
!! detovd = determinant of the overlap matrix where for the band
!!          that is being updated <C| is replaced by <D| (search direction)
!! dhc = Re[<D|H_0|C>]
!! dhd = <D|H_0|D>
!! efield_dot = reciprocal lattice coordinates of the electric field
!! hel(ifor,idir) = helicity of the ellipse associated with (ifor,idir)
!! nkpt = number of k-points
!! nsppol = 1 for unpolarized, 2 for spin-polarized
!! nstr(idir) = number of strings along the idir-th direction
!! sdeg = spin degeneracy
!! theta = value of the angle for which the energy (e0) and its
!!         derivative (e1) are computed
!!
!! OUTPUT
!! e0 = energy for the given value of theta
!! e1 = derivative of the energy with respect to theta
!!
!! PARENTS
!!      cgwf,linemin
!!
!! CHILDREN
!!      rhophi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine etheta(bcut,chc,detovc,detovd,dhc,dhd,efield_dot,e0,e1,&
&    hel,nkpt,nstr,sdeg,theta)

 use defs_basis
 use m_profiling_abi

 use m_numeric_tools, only : rhophi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'etheta'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt
 real(dp),intent(in) :: chc,dhc,dhd,sdeg,theta
 real(dp),intent(out) :: e0,e1
!arrays
 integer,intent(in) :: hel(2,3),nstr(3)
 real(dp),intent(in) :: bcut(2,3),detovc(2,2,3),detovd(2,2,3),efield_dot(3)

!Local variables -------------------------
!scalars
 integer :: idir,ifor
 real(dp) :: c2theta,ctheta,dphase,gnorm,phase,rho,s2theta,sgn,stheta
!arrays
 real(dp) :: dg_theta(2),g_theta(2)

! ***********************************************************************

 e0 = zero ; e1 = zero

 ctheta = cos(theta)
 stheta = sin(theta)
 c2theta = ctheta*ctheta - stheta*stheta   ! cos(2*theta)
 s2theta = two*ctheta*stheta               ! sin(2*theta)

 e0 = chc*ctheta*ctheta + dhd*stheta*stheta + dhc*s2theta
 e0 = e0*sdeg/nkpt

!DEBUG
!e0 = zero
!ENDDEBUG

 e1 = (dhd - chc)*s2theta + two*dhc*c2theta
 e1 = e1*sdeg/nkpt

 sgn = -1_dp
 do idir = 1, 3

   if (abs(efield_dot(idir)) < tol12) cycle

   do ifor = 1, 2

     g_theta(:)  = ctheta*detovc(:,ifor,idir) + &
&     stheta*detovd(:,ifor,idir)
     dg_theta(:) = -1_dp*stheta*detovc(:,ifor,idir) + &
&     ctheta*detovd(:,ifor,idir)

!    Compute E(theta)

     call rhophi(g_theta,phase,rho)
     if (theta >= bcut(ifor,idir)) phase = phase + hel(ifor,idir)*two_pi

!    DEBUG
!    unit = 100 + 10*idir + ifor
!    write(unit,'(4(f16.9))')theta,g_theta(:),phase
!    ENDDEBUG

     e0 = e0 + sgn*sdeg*efield_dot(idir)*phase/(two_pi*nstr(idir))


!    Compute dE/dtheta

!    imaginary part of the derivative of ln(g_theta)
     gnorm = g_theta(1)*g_theta(1) + g_theta(2)*g_theta(2)
     dphase = (dg_theta(2)*g_theta(1) - dg_theta(1)*g_theta(2))/gnorm

     e1 = e1 + sgn*sdeg*efield_dot(idir)*dphase/(two_pi*nstr(idir))

     sgn = -1_dp*sgn

   end do
 end do

end subroutine etheta
!!***

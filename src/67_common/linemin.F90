!{\src2tex{textfont=tt}}
!!****f* ABINIT/linemin
!! NAME
!! linemin
!!
!! FUNCTION
!! Performs the "line minimization" w.r.t. the angle theta on a unit circle
!! to update the wavefunction associated with the current k-point and
!! band label.
!! This routine is used only when the electric field is on (otherwise it could
!! in principle also be used, but there is a simpler procedure, as originally
!! coded in abinit).
!!
!! COPYRIGHT
!! Copyright (C) 2000-2017 ABINIT  group (MVeithen,ISouza,JIniguez)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! chc = <C|H_0|C> where |C> is the wavefunction of the current band
!! detovc = determinant of the overlap matrix S
!! detovd = determinant of the overlap matrix where for the band
!!          that is being updated <C| is replaced by <D| (search direction)
!! dhc = Re[<D|H_0|C>]
!! dhd = <D|H_0|D>
!! efield_dot = reciprocal lattice coordinates of the electric field
!! iline = index of the current line minimization
!! nkpt = number of k-points
!! nstr(idir) = number of strings along the idir-th direction
!! sdeg = spin degeneracy
!!
!! OUTPUT
!! bcut(ifor,idir) = branch cut of the ellipse associated with (ifor,idir)
!! costh = cos(thetam)
!! hel(ifor,idir) = helicity of the ellipse associated with (ifor,idir)
!! phase_end = total change in Zak phase, must be equal to
!!             dphase_aux1 + n*two_pi
!! sinth = sin(thetam)
!! thetam = optimal angle theta in line_minimization
!!
!!
!! SIDE EFFECTS
!! Input/Output
!! dphase_aux1 = change in Zak phase accumulated during the loop over iline
!!               (can be used for debugging in cgwf.f)
!! phase_init = initial Zak phase (before doing the first line minimization)
!!
!!
!! TODO
!!
!! NOTES
!! We are making the "frozen Hamiltonian approximation", i.e., the
!! Hamiltonian does not change with theta (we are neglecting the dependence
!! of the Hartree and exchange-correlation terms on theta; the original
!! abinit routine does the same)
!!
!! PARENTS
!!      cgwf
!!
!! CHILDREN
!!      etheta,rhophi,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine linemin(bcut,chc,costh,detovc,detovd,dhc,dhd,dphase_aux1,&
&  efield_dot,iline,nkpt,nstr,hel,phase_end,phase_init,sdeg,sinth,thetam)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_numeric_tools, only : rhophi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'linemin'
 use interfaces_14_hidewrite
 use interfaces_67_common, except_this_one => linemin
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iline,nkpt
 real(dp),intent(in) :: chc,dhc,dhd,sdeg
 real(dp),intent(out) :: costh,sinth,thetam
!arrays
 integer,intent(in) :: nstr(3)
 integer,intent(out) :: hel(2,3)
 real(dp),intent(in) :: detovc(2,2,3),detovd(2,2,3),efield_dot(3)
 real(dp),intent(inout) :: dphase_aux1(3),phase_init(3)
 real(dp),intent(out) :: bcut(2,3),phase_end(3)

!Local variables -------------------------
!scalars
 integer :: idir,ifor,igrid,iter,maxiter,ngrid
 real(dp) :: aa,angle,bb,big_axis,cc,cphi_0,delta_theta,e0,e1
 real(dp) :: excentr,iab,phase0,phase_min,phi_0,rdum,sgn,small_axis,sphi_0
 real(dp) :: theta,theta_0,val
 logical :: flag_neg
 character(len=500) :: message
!arrays
 real(dp) :: g_theta(2),theta_min(2),theta_try(2)
 real(dp) :: esave(251),e1save(251)   !!REC

! ***********************************************************************

!Compute the helicity and the branch cut of the ellipse in the complex
!plane associated with the overlap between a k-point and one of its neighbours

 do idir = 1, 3

   if (abs(efield_dot(idir)) < tol12) cycle

   do ifor = 1, 2

     aa = half*(detovc(1,ifor,idir)*detovc(1,ifor,idir) + &
&     detovc(2,ifor,idir)*detovc(2,ifor,idir) + &
&     detovd(1,ifor,idir)*detovd(1,ifor,idir) + &
&     detovd(2,ifor,idir)*detovd(2,ifor,idir))

     bb = half*(detovc(1,ifor,idir)*detovc(1,ifor,idir) + &
&     detovc(2,ifor,idir)*detovc(2,ifor,idir) - &
&     detovd(1,ifor,idir)*detovd(1,ifor,idir) - &
&     detovd(2,ifor,idir)*detovd(2,ifor,idir))

     cc = detovc(1,ifor,idir)*detovd(1,ifor,idir) + &
&     detovc(2,ifor,idir)*detovd(2,ifor,idir)

     iab = detovc(1,ifor,idir)*detovd(2,ifor,idir) - &
&     detovc(2,ifor,idir)*detovd(1,ifor,idir)

     if (iab >= zero) then
       hel(ifor,idir) = 1
     else
       hel(ifor,idir) = -1
     end if

     if (abs(bb) > tol8) then
       theta_0 = half*atan(cc/bb)
     else
       theta_0 = quarter*pi
     end if

     if (bb < zero) theta_0 = theta_0 + pi*half

     g_theta(:) = cos(theta_0)*detovc(:,ifor,idir) + &
&     sin(theta_0)*detovd(:,ifor,idir)
!    DEBUG
!    write(std_out,*)'before rhophi, g_theta =',g_theta
!    ENDDEBUG
     call rhophi(g_theta,phi_0,rdum)
!    DEBUG
!    write(std_out,*)'after rhophi, phi_0 = ',phi_0
!    ENDDEBUG

     cphi_0 = cos(phi_0)
     sphi_0 = sin(phi_0)

     rdum = aa - sqrt(bb*bb + cc*cc)
     if (rdum < zero) rdum = zero
     small_axis = sqrt(rdum)
     big_axis = sqrt(aa + sqrt(bb*bb + cc*cc))
     excentr = hel(ifor,idir)*small_axis/big_axis

!    Find angle for which phi = pi
     if (abs(excentr) > tol8) then
       angle = atan(tan(pi-phi_0)/excentr)
     else
       if (tan(pi-phi_0)*hel(ifor,idir) > zero) then
         angle = half*pi
       else
         angle = -0.5_dp*pi
       end if
     end if
     bcut(ifor,idir) = angle + theta_0


!    Compute the branch-cut angle
     if (hel(ifor,idir) == 1) then
       if ((sphi_0 > 0).and.(cphi_0 > 0)) bcut(ifor,idir) = bcut(ifor,idir) + pi
       if ((sphi_0 < 0).and.(cphi_0 > 0)) bcut(ifor,idir) = bcut(ifor,idir) - pi
     else
       if ((sphi_0 > 0).and.(cphi_0 > 0)) bcut(ifor,idir) = bcut(ifor,idir) - pi
       if ((sphi_0 < 0).and.(cphi_0 > 0)) bcut(ifor,idir) = bcut(ifor,idir) + pi
     end if

     if (bcut(ifor,idir) > pi) bcut(ifor,idir) = bcut(ifor,idir) - two_pi
     if (bcut(ifor,idir) < -1_dp*pi) bcut(ifor,idir) = bcut(ifor,idir) + two_pi

!    DEBUG
!    write(std_out,'(a,2x,i3,2x,i3,5x,f16.9,5x,i2)')'linemin: ifor,idir,bcut,hel',&
!    &   ifor,idir,bcut(ifor,idir),hel(ifor,idir)
!    write(std_out,'(a,5x,f16.9,5x,f16.9)')'linemin: big_axis,small_axis ',&
!    &     big_axis,small_axis
!    ENDDEBUG

   end do   ! ifor
 end do   ! idir

!---------------------------------------------------------------------------

!Perform the "line minimization" w.r.t. the angle theta on a unit circle
!to update the wavefunction associated with the current k-point and band label.

 ngrid = 250   ! initial number of subdivisions in [-pi/2,pi/2]
!for finding extrema
 maxiter = 100
 delta_theta = pi/ngrid

!DEBUG
!write(std_out,*)'linemin: theta, e0, e1, e1fdiff'
!ENDDEBUG


!Get the interval where the absolute minimum of E(theta) is located

 val = huge(one)             ! large number
 flag_neg=.false.
 theta_min(:) = ten
 do igrid = 1, ngrid+1

   theta = (igrid - 1)*delta_theta - pi*half
   call etheta(bcut,chc,detovc,detovd,dhc,dhd,efield_dot,e0,e1,&
&   hel,nkpt,nstr,sdeg,theta)

   esave(igrid)=e0      !!REC
   e1save(igrid)=e1     !!REC

!  It is important to detect when the slope changes from negative to positive
!  Moreover, a slope being extremely close to zero must be ignored

!  DEBUG
!  write(std_out,*)' igrid,e0,e1,val,theta_min(:)=',igrid,theta,e0,e1,val,theta_min(:)
!  ENDDEBUG
   
!  Store e1 and theta if negative ...
   if(e1 < -tol10)then
     theta_try(1)=theta
     flag_neg=.true.
   end if
!  A change of sign is just happening
   if(e1 > tol10 .and. flag_neg)then
     theta_try(2)=theta
     flag_neg=.false.
!    Still, must be better than the previous minimum in order to succeed
     if (e0 < val-tol10) then
       val=e0
       theta_min(:)=theta_try(:)
     end if
   end if
 end do

!In case the minimum was not found

 if (abs(theta_min(1) - ten) < tol10) then
!  REC start
   write(message,'(a,a)')ch10,&
&   ' linemin: ERROR- cannot find theta_min.'
   call wrtout(std_out,message,'COLL')
   write(message,'(a,a)')ch10,&
&   ' igrid      theta          esave(igrid)    e1save(igrid) '
   call wrtout(std_out,message,'COLL')
   do igrid = 1, ngrid+1
     theta = (igrid - 1)*delta_theta - pi*half
     write(std_out,'(i6,3f16.9)')igrid,theta,esave(igrid),e1save(igrid)
     !write(101,'(i6,3f16.9)')igrid,theta,esave(igrid),e1save(igrid)
   end do
   write(message,'(6a)')ch10,&
&   ' linemin : ERROR - ',ch10,&
&   '  Cannot find theta_min. No minimum exists : the field is too strong ! ',ch10,&
&   '  Try decreasing difference between D and 4 Pi P by changing structure or D (only for fixed D calculation)'
   call wrtout(std_out,message,'COLL')
   message = ' linemin cannot find theta_min'
   MSG_ERROR(message)
 end if

!Compute the mimum of E(theta)


 iter = 0
 do while ((delta_theta > tol8).and.(iter < maxiter))
   delta_theta = half*(theta_min(2) - theta_min(1))
   theta = theta_min(1) + delta_theta
   call etheta(bcut,chc,detovc,detovd,dhc,dhd,efield_dot,e0,e1,&
&   hel,nkpt,nstr,sdeg,theta)
   if (e1 > zero) then
     theta_min(2) = theta
   else
     theta_min(1) = theta
   end if
   iter = iter + 1

!  DEBUG
!  write(std_out,'(a,2x,i3,2(2x,f16.9))')'iter,e0,e1 = ',iter,e0,e1
!  ENDDEBUG

 end do

 costh = cos(theta)
 sinth = sin(theta)

 thetam = theta

!DEBUG
!write(std_out,*)'linemin : thetam = ',thetam
!ENDDEBUG


!---------------------------------------------------------------------------

!Compute and store the change in electronic polarization

 sgn = one
 do idir = 1, 3

   if (abs(efield_dot(idir)) < tol12) cycle

   phase_end(idir) = zero
   do ifor = 1, 2

     g_theta(:) = detovc(:,ifor,idir)
!    DEBUG
!    write(std_out,*)'before rhophi (2nd call), g_theta =',g_theta
!    ENDDEBUG
     call rhophi(g_theta,phase0,rdum)
!    DEBUG
!    write(std_out,*)'after rhophi, phase0 = ',phase0
!    ENDDEBUG

     if(iline == 1) phase_init(idir) = phase_init(idir) + sgn*phase0

     g_theta(:) = costh*detovc(:,ifor,idir) + sinth*detovd(:,ifor,idir)
     call rhophi(g_theta,phase_min,rdum)

     phase_end(idir) = phase_end(idir) + sgn*phase_min

!    Correct for branch cuts (remove jumps)
     if (bcut(ifor,idir) <= zero) phase0 = phase0 + hel(ifor,idir)*two_pi
     if(thetam >= bcut(ifor,idir)) phase_min = phase_min + hel(ifor,idir)*two_pi

     dphase_aux1(idir) = dphase_aux1(idir) + sgn*(phase_min - phase0)

     sgn = -1_dp*sgn

   end do   ! idir
 end do    ! ifor

!DEBUG
!write(std_out,'(a,3(2x,f16.9))')'dphase_aux1 = ',(dphase_aux1(idir),idir = 1, 3)
!write(std_out,*)' linemin: debug, exit.'
!ENDDEBUG

end subroutine linemin
!!***

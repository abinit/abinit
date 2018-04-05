!{\src2tex{textfont=tt}}
!!****f* ABINIT/rsurf
!! NAME
!! rsurf
!!
!! FUNCTION
!! Basic routine for determination of the radius of Bader surface
!! for spherical rayon theta,phi
!! the bassin is tested by following the gradient line
!! If srch==true (in general for calls from surf) the routine aim_follow
!! is called to stop when it arrives under already known part of surface
!! Simple bissection method is used to obtain the radius
!!
!! COPYRIGHT
!! Copyright (C) 2002-2018 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! INPUTS
!! aim_dtset= the structured entity containing all input variables
!! rr0= starting radius
!! theta,phi = the spherical direction
!! iatinit= the atom index
!! srch= see above
!! npmax= maximum number of divisions in one step for follow
!!
!! OUTPUT
!! rr= radius
!! grho(3)= gradient on the surface
!!
!! PARENTS
!!      drvaim,surf
!!
!! CHILDREN
!!      aim_follow,timein,vgh_rho
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine rsurf(aim_dtset,rr,grho,theta,phi,rr0,iatinit,npmax,srch)

 use m_profiling_abi

 use defs_basis
 use defs_parameters
 use defs_aimprom
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rsurf'
 use interfaces_18_timing
 use interfaces_63_bader, except_this_one => rsurf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iatinit,npmax
 real(dp),intent(in) :: phi,rr0,theta
 real(dp),intent(out) :: rr
 logical,intent(in) :: srch
!arrays
 real(dp),intent(out) :: grho(3)
!no_abirules
 type(aim_dataset_type),intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: iat,ii,ipos,iposinit,jj,nstep
 real(dp),parameter :: mfkt=1.d1
 real(dp) :: aa,dmax,dr,drr,rho,rr1,rr2,t1,t2,wall
 logical :: cross,deb_tmp,in,in1,in2,low,srch_tmp
!arrays
 real(dp) :: hrho(3,3),unvec(3),vv(3)

! *********************************************************************

 srch_tmp=srch
 deb_tmp=deb

!unity vecteur in the direction (theta,phi)

 unvec(1)=sin(theta)*cos(phi)
 unvec(2)=sin(theta)*sin(phi)
 unvec(3)=cos(theta)


 rr=rr0
 rr1=rr
 rr2=rr
 drr=1._dp
 if (abs(rr0-r0)<1.0d-12) then
   dr=aim_dtset%dr0*mfkt
 else
   dr=aim_dtset%dr0
 end if

 vv(1)=xatm(1,aim_dtset%batom)
 vv(2)=xatm(2,aim_dtset%batom)
 vv(3)=xatm(3,aim_dtset%batom)


 iposinit=batcell
 write(std_out,'("ATOM iat=",i4," ipos=",i4)') aim_dtset%batom,batcell
 jj=0

 cross=.false.

 in=.true.
 low=.false.

 dmax=h0

 in1=.true.
 in2=in1

 do while((drr>aim_drmin).or.(jj<2))
   call timein(t1,wall)
   jj=jj+1
   do ii=1,3
     vv(ii)=xatm(ii,aim_dtset%batom)+rr*unvec(ii)
   end do

!  VACUUM CONDITION

   call vgh_rho(vv,rho,grho,hrho,aa,iat,ipos,0)
   if (rho < aim_rhomin) exit

   ldeb=.false.

   call aim_follow(aim_dtset,vv,npmax,srch_tmp,iatinit,iposinit,iat,ipos,nstep)

   call timein(t2,wall)
   t2=t2-t1

   write(std_out,'(a,i4,a,f12.8,a,i4,a,i4,a,f10.5,a,i4)') &
&   ' :STEP ',jj,' r=',rr,' iat=',iat,' ipos=',ipos,' time(sec)=',t2,' nstep=',nstep

   if ((iat.eq.iatinit).and.(ipos.eq.iposinit)) then
     in=.true.
   else
     in=.false.
   end if

!  
!  NEW RADIUS
!  

   if ((jj.eq.1).or.((in1.eqv.in).and.(.not.cross))) then
     if (in) then
       rr2=rr1
       rr1=rr
       rr=rr+dr
     else
       rr2=rr1
       rr1=rr
       rr=rr-dr
     end if
     if ((jj>2).and.(dr<(0.6))) then
!      modification of the step
       dr=dr*aim_fac
       if (deb_tmp) write(std_out,*) ':DR ',dr
     end if
   else
     if (.not.cross) then
       cross=.true.
       rr2=rr1
     else
       if (in2) then
         if (in) then
           rr2=rr1
         else
           in1=in2
         end if
       else
         if (in) then
           in1=in2
         else
           rr2=rr1
         end if
       end if
     end if
     rr1=rr
     rr=(rr2+rr1)/2.0
   end if

   in2=in1
   in1=in
   drr=abs(rr2-rr1)/rr
   if (deb_tmp) write(std_out,*) ':DRR ',jj,rr2,rr1,drr
 end do

end subroutine rsurf
!!***

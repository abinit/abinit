!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_sphharm
!! NAME
!!  m_paw_sphharm
!!
!! FUNCTION
!!  This module contains a set of routines to compute the complex (resp. real)
!!  spherical harmonics Ylm (resp. Slm) (and gradients).
!!
!! COPYRIGHT
!! Copyright (C) 2013-2017 ABINIT group (MT, FJ, TRangel)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

MODULE m_paw_sphharm

 USE_DEFS
 USE_MSG_HANDLING
 USE_MEMORY_PROFILING

 implicit none

 private

!public procedures.
 public :: ylmc         ! Complex Spherical harmonics for l<=3.
 public :: ylmcd        ! First derivative of complex Ylm wrt theta and phi up to l<=3
 public :: ylm_cmplx    ! All (complex) spherical harmonics for lx<=4
 public :: initylmr     ! Real Spherical Harmonics on a set of vectors
 public :: ys           ! Matrix element <Yl'm'|Slm>
 public :: lxyz         ! Matrix element <Yl'm'|L_idir|Ylm>
 public :: slxyzs       ! Matrix element <Sl'm'|L_idir|Slm>
 public :: plm_coeff    ! Coefficients depending on Plm used to compute the 2nd der of Ylm
 public :: ass_leg_pol  ! Associated Legendre Polynomial Plm(x)
 public :: plm_d2theta  ! d2(Plm (cos(theta)))/d(theta)2 (P_lm= ass. legendre polynomial)
 public :: plm_dphi     ! m*P_lm(x)/sqrt((1-x^2)  (P_lm= ass. legendre polynomial)
 public :: plm_dtheta   ! -(1-x^2)^1/2*d/dx{P_lm(x)} (P_lm= ass. legendre polynomial)
 public :: pl_deriv     ! d2(Pl (x)))/d(x)2  where P_l is a legendre polynomial
 public :: mat_mlms2jmj ! Change a matrix from the Ylm basis to the J,M_J basis
 public :: mat_slm2ylm  ! Change a matrix from the Slm to the Ylm basis or from Ylm to Slm
!!***

CONTAINS

!===========================================================
!!***

!!****f* m_paw_sphharm/ylmc
!! NAME
!! ylmc
!!
!! FUNCTION
!!  Return a complex spherical harmonic with l <= 3
!!
!! INPUTS
!!  il=angular quantum number
!!  im=magnetic quantum number
!!  kcart=vector in cartesian coordinates defining the value of \theta and \psi
!!   where calculate the spherical harmonic
!!
!! OUTPUT
!!  ylm= spherical harmonic
!!
!! NOTES
!!  Note the use of double precision complex.
!!  Case l>3 not implemented.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ylmc(il,im,kcart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ylmc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: il,im
 complex(dpc) :: ylmc
!arrays
 real(dp),intent(in) :: kcart(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: LMAX=3
 real(dp),parameter :: PPAD=tol8
 real(dp) :: cosphi,costh,costhreephi,costwophi,r,rxy,sinphi,sinth,sinthreephi,sintwophi
 !complex(dpc) :: new_ylmc
 character(len=500) :: msg

! *************************************************************************

 if (ABS(im)>ABS(il)) then
   write(msg,'(3(a,i0))') 'm is,',im,' however it should be between ',-il,' and ',il
   MSG_ERROR(msg)
 end if

 ylmc = czero

 r=SQRT(kcart(1)**2+kcart(2)**2+kcart(3)**2)
 if (r<PPAD) r=r+PPAD
 !$if (r<tol10) RETURN

 rxy=SQRT(kcart(1)**2+kcart(2)**2)
 if (rxy<PPAD)rxy=r+PPAD
!
! Determine theta and phi 
 costh= kcart(3)/r

#if 1
 ! old buggy coding
 sinth= rxy/r
 cosphi= kcart(1)/rxy
 sinphi= kcart(2)/rxy
#else
 sinth=sqrt(abs((one-costh)*(one+costh))) ! abs is needed to prevent very small negative arg
 cosphi=one
 sinphi=zero
 if (sinth>tol10) then
   cosphi=kcart(1)/(r*sinth)
   sinphi=kcart(2)/(r*sinth)
 end if
#endif

 costwophi= two*cosphi**2 - one
 sintwophi= two*sinphi*cosphi
 costhreephi=cosphi*costwophi-sinphi*sintwophi
 sinthreephi=cosphi*sintwophi+sinphi*costwophi

 select case (il)

 case (0)
  ylmc= one/SQRT(four_pi)

 case (1)
  if (ABS(im)==0) then
   ylmc = SQRT(three/(four_pi))*costh
  else if (ABS(im)==1) then
   ylmc = -SQRT(three/(8._dp*pi))*sinth*CMPLX(cosphi,sinphi)
  else
   msg='wrong im'
   MSG_ERROR(msg)
  end if

 case (2)
  if (ABS(im)==0) then
   ylmc = SQRT(5.d0/(16.d0*pi))*(three*costh**2-one)
  else if (ABS(im)==1) then
   ylmc = -SQRT(15.d0/(8.d0*pi))*sinth*costh*cmplx(cosphi,sinphi)
  else if (ABS(im)==2) then
   ylmc = SQRT(15.d0/(32.d0*pi))*(sinth)**2*CMPLX(costwophi,sintwophi)
  else 
   msg='wrong im'
   MSG_ERROR(msg)
  end if

 case (3)
  if (ABS(im)==0) then
   ylmc= SQRT(7.d0/(16.d0*pi))*(5.d0*costh**3 -3.d0*costh)
  else if (ABS(im)==1) then
   ylmc= -SQRT(21.d0/(64.d0*pi))*sinth*(5.d0*costh**2-one)*CMPLX(cosphi,sinphi)
  else if (ABS(im)==2) then
   ylmc= SQRT(105.d0/(32.d0*pi))*sinth**2*costh*CMPLX(costwophi,sintwophi)
  else if (ABS(im)==3) then
   ylmc=-SQRT(35.d0/(64.d0*pi))*sinth**3*CMPLX(costhreephi,sinthreephi)
  else 
   msg='wrong im'
   MSG_ERROR(msg)
  end if

 case default
  !write(msg,'(a,i6,a,i6)')' The maximum allowed value for l is,',LMAX,' however l=',il
  !MSG_ERROR(msg)
 end select
!
!=== Treat the case im < 0 ===
 if (im < 0) ylmc=(-one)**(im)*CONJG(ylmc)

 ! FIXME: Use the piece of code below as it works for arbitrary (l,m)
 ! the implementation above is buggy when the vector is along z! 
 ! 
#if 0
! Remember the expression of complex spherical harmonics:
! $Y_{lm}(\theta,\phi)=sqrt{{(2l+1) over (4\pi)} {fact(l-m)/fact(l+m)} } P_l^m(cos(\theta)) e^{i m\phi}$
  new_ylmc = SQRT((2*il+1)*factorial(il-ABS(im))/(factorial(il+ABS(im))*four_pi)) * &
&   ass_leg_pol(il,ABS(im),costh) * CMPLX(cosphi,sinphi)**ABS(im)
  if (im<0) new_ylmc=(-one)**(im)*CONJG(new_ylmc)

  if (ABS(new_ylmc-ylmc)>tol6) then
    !MSG_WARNING("Check new_ylmc")
    !write(std_out,*)"il,im,new_ylmc, ylmc",il,im,new_ylmc,ylmc
    !write(std_out,*)"fact",SQRT((2*il+1)*factorial(il-ABS(im))/(factorial(il+ABS(im))*four_pi))
    !write(std_out,*)"costh,sinth,ass_leg_pol",costh,sinth,ass_leg_pol(il,ABS(im),costh) 
    !write(std_out,*)"cosphi,sinphi,e^{imphi}",cosphi,sinphi,CMPLX(cosphi,sinphi)**ABS(im)
  end if
  ylmc = new_ylmc
#endif

end function ylmc
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/ylmcd
!! NAME
!! ylmcd
!!
!! FUNCTION
!!  Computes dth and dphi, the first derivatives of complex Ylm as a function of 
!!  th and phi (the angles of the spherical coordinates)
!!  It works for all spherical harmonics with l <= 3
!!
!! INPUTS
!!  il=angular quantum number
!!  im=magnetic quantum number
!!  kcart=cartesian coordinates of the vector where the first derivatives of Ylm are evaluated
!!
!! OUTPUT
!!  dth =derivative of Y_lm with respect to \theta
!!  dphi=derivative of Y_lm with respect to \phi
!!
!! NOTES
!!  Note the use of double precision complex.
!!  Case l>3 not implemented.
!!
!! PARENTS
!!      m_vkbr
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine ylmcd(il,im,kcart,dth,dphi)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ylmcd'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: il,im
 complex(dpc),intent(out) :: dphi,dth
!arrays
 real(dp),intent(in) :: kcart(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: LMAX=3
 real(dp),parameter :: PPAD=tol8
 real(dp) :: cosphi,costh,costhreephi,costwophi,r,rxy,sinphi,sinth,sinthreephi,sintwophi
 character(len=500) :: msg

! *************************************************************************

 if (ABS(im)>ABS(il))then
   write(msg,'(3(a,i0))')' m is,',im,' however it should be between ',-il,' and ',il
   MSG_ERROR(msg)
 end if

 dphi=czero; dth=czero

 r=SQRT(kcart(1)**2+kcart(2)**2+kcart(3)**2)
 if (r<PPAD) r=r+PPAD
 !$if (r<tol10) RETURN

 rxy=SQRT(kcart(1)**2+kcart(2)**2)
 if (rxy<PPAD) rxy=r+PPAD

! Determine theta and phi 
 costh= kcart(3)/r
#if 1
 ! old buggy coding
 sinth= rxy/r
 cosphi= kcart(1)/rxy
 sinphi= kcart(2)/rxy
#else
 sinth=sqrt(abs((one-costh)*(one+costh))) ! abs is needed to prevent very small negative arg
 cosphi=one
 sinphi=zero
 if (sinth>tol10) then
   cosphi=kcart(1)/(r*sinth)
   sinphi=kcart(2)/(r*sinth)
 end if
#endif

 costwophi= two*cosphi**2 - one
 sintwophi= two*sinphi*cosphi
 costhreephi=cosphi*costwophi-sinphi*sintwophi
 sinthreephi=cosphi*sintwophi+sinphi*costwophi

 select case (il)

 case (0)
   dth  = czero 
   dphi = czero

 case (1)
   if (ABS(im)==0) then
     dth= -SQRT(three/(four_pi))*sinth
     dphi= czero
   else if (abs(im)==1) then
     dth= -SQRT(3.d0/(8.d0*pi))*costh*CMPLX(cosphi,sinphi)
     dphi=-SQRT(3.d0/(8.d0*pi))*sinth*CMPLX(-sinphi,cosphi)
   end if

 case (2)
   if (ABS(im)==0) then
     dth= -SQRT(5.d0/(16.d0*pi))*6.d0*costh*sinth
     dphi= czero
   else if (ABS(im)==1) then
     dth=  -SQRT(15.d0/(8.d0*pi))*(costh**2-sinth**2)*CMPLX(cosphi,sinphi)
     dphi= -SQRT(15.d0/(8.d0*pi))*costh*sinth*(0.d0,1.d0)*CMPLX(cosphi,sinphi)
   else if (abs(im)==2) then
     dth  = SQRT(15.d0/(32.d0*pi))*2.d0*costh*sinth*CMPLX(costwophi,sintwophi)
     dphi = SQRT(15.d0/(32.d0*pi))*sinth**2*(0.d0,2.d0)*CMPLX(costwophi,sintwophi)
   end if

 case (3)
   if (ABS(im)==0) then
     dth = SQRT(7.d0/(16*pi))*(-15.d0*costh**2*sinth + 3.d0**sinth)
     dphi= czero
   else if (ABS(im)==1) then
     dth= -SQRT(21.d0/(64.d0*pi))*CMPLX(cosphi,sinphi)*(5.d0*costh**3-costh-10.d0*sinth**2*costh)
     dphi=-SQRT(21.d0/(64.d0*pi))*sinth*(5.d0*costh**2-1)*(0.d0,1.d0)*CMPLX(cosphi,sinphi)
   else if (ABS(im)==2) then
     dth =SQRT(105.d0/(32.d0*pi))*(2.d0*sinth*costh**2-sinth**3)*CMPLX(costwophi,sintwophi)
     dphi=SQRT(105.d0/(32*pi))*sinth**2*costh*(0.d0,2.d0)*CMPLX(costwophi,sintwophi)
   else if (abs(im)==3) then
     dth =-SQRT(35.d0/(64.d0*pi))*3.d0*sinth**2*costh*CMPLX(costhreephi,sinthreephi)
     dphi= SQRT(35.d0/(64.d0*pi))*sinth**3*(0.d0,3.d0)*CMPLX(costhreephi,sinthreephi)
   end if

 case default
   write(msg,'(2(a,i0))')' The maximum allowed value for l is,',LMAX,' however, l=',il
   MSG_ERROR(msg)
 end select
!
!=== Treat the case im < 0 ===
 if (im<0) then
   dth = (-one)**(im)*CONJG(dth)
   dphi= (-one)**(im)*CONJG(dphi)
 end if

end subroutine ylmcd
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/ylm_cmplx
!! NAME
!! ylm_cmplx
!!
!! FUNCTION
!!  Calculate all (complex) spherical harmonics for lx<=4
!!
!! INPUTS
!!  lx= quantum numbers.
!!  xx= cartesian coordinate in the x direction
!!  yy= cartesian coordinate in the y direction
!!  zz= cartesian coordinate in the z direction
!!
!! cartesian coordinates
!! OUTPUT
!!  ylm((lx+1)*(lx+1)) complex spherical harmonics for all l<=lx and all
!!                     possible values of m.
!!
!! NOTES
!!  We are supressing the so-called Condon-Shortley phase
!!
!! PARENTS
!!      mlwfovlp_proj,mlwfovlp_ylmfac
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine ylm_cmplx(lx,ylm,xx,yy,zz)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ylm_cmplx'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lx
 real(dp),intent(in) :: xx,yy,zz
!arrays
 complex(dpc),intent(out) :: ylm((lx+1)*(lx+1))

!Local variables-------------------------------
!scalars
 integer :: ii,l1,m1,nc,nn
 real(dp) :: dc,dl,dm,ds,rr,rrs,rs,sq2,w,x,xs,ya,yi,yr
!arrays
 real(dp) :: cosa(lx+1),fact(2*(lx+1)),plm(lx+2,lx+2),qlm(lx+2,lx+2),sgn(lx+1)
 real(dp) :: sina(lx+1)

! *************************************************************************
 
!normalization coefficients
 sq2=sqrt(2.0d0)
 fact(1)=1.0d0
 do ii=2,2*lx+1
   fact(ii)=(ii-1)*fact(ii-1)
 end do
 do l1=1,lx+1
   sgn(l1)=(-1.d0)**(l1-1)
   do m1=1,l1
     qlm(l1,m1)=sqrt((2*l1-1)*fact(l1-m1+1)/&
&     (four_pi*fact(l1+m1-1)))
   end do
 end do

!legendre polynomials
 rs=xx**2 + yy**2 + zz**2
 if(rs > tol8) then
   xs=zz**2/rs
   x=zz/sqrt(rs)
   w=sqrt(abs(1.0d0 - xs))
 else
   x=0.0d0

   w=1.0d0
 end if
 plm(1,1)=1.0d0
 plm(2,1)=x
 plm(2,2)=w
 plm(3,2)=3.0d0*x*w
 do m1=1,lx
   dm=m1-1
   if(m1 > 1) then
     plm(m1+1,m1)=x*plm(m1,m1) + 2*dm*w*plm(m1,m1-1)
   end if
   if(m1 < lx) then
     do l1=m1+2,lx+1
       dl=l1-1
       plm(l1,m1)=((2*dl-1)*x*plm(l1-1,m1)&
&       - (dl+dm-1)*plm(l1-2,m1))/(dl-dm)
     end do
   end if
   plm(m1+1,m1+1)=(2*dm+1)*w*plm(m1,m1)
 end do

!azimuthal angle phase factors
 rrs=xx**2 + yy**2
 if(rrs > tol8) then
   rr=sqrt(rrs)
   dc=xx/rr
   ds=yy/rr
 else
   dc=1.0d0
   ds=0.0d0
 end if
 cosa(1)=1.0d0
 sina(1)=0.0d0
 do m1=2,lx+1
   cosa(m1)=dc*cosa(m1-1) - ds*sina(m1-1)
   sina(m1)=ds*cosa(m1-1) + dc*sina(m1-1)
 end do

!combine factors
 do l1=1,lx+1
   do m1=2,l1
     nn=(l1-1)**2 + (l1-1) + (m1-1) + 1
     nc=(l1-1)**2 + (l1-1) - (m1-1) + 1
!    note that we are supressing the so-called Condon-Shortley phase
!    ya=sgn(m1)*qlm(l1,m1)*plm(l1,m1)
     ya=qlm(l1,m1)*plm(l1,m1)
     yr=ya*cosa(m1)
     yi=ya*sina(m1)
     ylm(nc)=sgn(m1)*cmplx(yr,-yi)
     ylm(nn)=cmplx(yr,yi)
   end do
 end do
 do l1=1,lx+1
   nn=(l1-1)**2 + (l1-1) + 1
   ya=qlm(l1,1)*plm(l1,1)
   ylm(nn)=cmplx(ya,0.d0)
 end do

end subroutine ylm_cmplx
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/initylmr
!! NAME
!! initylmr
!!
!! FUNCTION
!! Calculate the real spherical harmonics Ylm (and gradients)
!! over a set of (r) vectors given in Cartesian coordinates.
!!
!! INPUTS
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotential
!!  normchoice=0  the input rr vectors are normalized
!!            =1  the norm of the input vector is in nrm() array
!!  nrm(npts) = Depending of normchoice, this array contains
!!              either the weight of the point or the norm of rr.
!!  npts = number of rr vectors
!!  option= 1=compute Ylm(R), 2=compute Ylm(R) and dYlm/dRi (cartesian derivatives),
!!          3=compute Ylm(R), dYlm/dRi and d2Ylm/dRidRj (cartesian derivatives)
!!  rr(3,npts)=  vectors for which ylmr have to be calculated
!!               For each point of the spherical mesh, gives the
!!               Cartesian coordinates of the corresponding point.
!!
!! OUTPUT
!!  if (option=1, 2 or 3)
!!    ylm(mpsang*mpsang,npts)     = real spherical harmonics for each r point
!!  if (option=2 or 3)
!!    ylmr_gr(1:3,mpsang*mpsang,npts)= gradients of real spherical harmonics
!!  if (option=3)
!!    ylmr_gr(4:9,mpsang*mpsang,npts)= first and second gradients of real spherical harmonics
!!
!! NOTES
!! Remember the expression of complex spherical harmonics:
!! $Y_{lm}(%theta ,%phi)=sqrt{{(2l+1) over (4 %pi)} {fact(l-m) over fact(l+m)} } P_l^m(cos(%theta)) func e^{i m %phi}$
!! Remember the expression of real spherical harmonics as linear combination of imaginary spherical harmonics:
!! $Yr_{lm}(%theta ,%phi)=(Re{Y_{l-m}}+(-1)^m Re{Y_{lm}})/sqrt{2}
!! $Yr_{l-m}(%theta ,%phi)=(Im{Y_{l-m}}-(-1)^m Im{Y_{lm}})/sqrt{2}
!!
!! PARENTS
!!      debug_tools,denfgr,m_paw_finegrid,m_paw_pwaves_lmn,m_pawang
!!      mlwfovlp_ylmfar,posdoppler,pspnl_operat_rec,qijb_kk,smatrix_pawinit
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine initylmr(mpsang,normchoice,npts,nrm,option,rr,ylmr,ylmr_gr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initylmr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpsang,normchoice,npts,option
!arrays
 real(dp),intent(in) :: nrm(npts),rr(3,npts)
 real(dp),intent(out) :: ylmr(mpsang*mpsang,npts)
 real(dp),optional,intent(out) :: ylmr_gr(3*(option/2)+6*(option/3),mpsang*mpsang,npts)

!Local variables ------------------------------
!scalars
 integer :: dimgr,ilang,inpt,l0,ll,mm
 real(dp) :: cphi,ctheta,fact,onem,rnorm,sphi,stheta,work1,work2,ylmcst,ylmcst2
 logical :: compute_ylm,compute_ylm2gr,compute_ylmgr
!arrays
 real(dp) :: dphi(3),dtheta(3),iphase(mpsang-1),rphase(mpsang-1)
 real(dp),allocatable :: blm(:,:)

!************************************************************************

!What has to be computed ?
 compute_ylm   = (option==1.or.option==2.or.option==3)
 compute_ylmgr =((             option==2.or.option==3).and.present(ylmr_gr))
 compute_ylm2gr=((                          option==3).and.present(ylmr_gr))
 dimgr=3*(option/2)+6*(option/3)

!Initialisation of spherical harmonics
 if (compute_ylm  ) ylmr   (:  ,1:npts)=zero
 if (compute_ylmgr) ylmr_gr(:,:,1:npts)=zero

!Special case for l=0
 if (compute_ylm  ) ylmr(1,1:npts)=1._dp/sqrt(four_pi)
 if (compute_ylmgr) ylmr_gr(1:dimgr,1,1:npts)=zero
 if (mpsang>1) then

!  Loop over all rr
   do inpt=1,npts

!    Load module of rr
     rnorm=one
     if (normchoice==1) rnorm=nrm(inpt)

!    Continue only for r<>0

     if (rnorm>tol10) then

!      Determine theta and phi
       cphi=one
       sphi=zero
       ctheta=rr(3,inpt)/rnorm
!      MM030519 : abs is needed to prevent very small negative arg
       stheta=sqrt(abs((one-ctheta)*(one+ctheta)))
       if (stheta>tol10) then
         cphi=rr(1,inpt)/(rnorm*stheta)
         sphi=rr(2,inpt)/(rnorm*stheta)
       end if
       do mm=1,mpsang-1
         rphase(mm)=dreal(dcmplx(cphi,sphi)**mm)
         iphase(mm)=aimag(dcmplx(cphi,sphi)**mm)
       end do

!      Determine gradients of theta and phi
       if (compute_ylmgr) then
         dtheta(1)=ctheta*cphi
         dtheta(2)=ctheta*sphi
         dtheta(3)=-stheta
         dphi(1)=-sphi
         dphi(2)=cphi
         dphi(3)=zero
       end if

!      COMPUTE Ylm(R)
       if (compute_ylm) then
!        Loop over angular momentum l
         do ilang=2,mpsang
           ll=ilang-1
           l0=ll**2+ll+1
           fact=1._dp/real(ll*(ll+1),dp)
           ylmcst=sqrt(real(2*ll+1,dp)/four_pi)
!          Special case m=0
           ylmr(l0,inpt)=ylmcst*ass_leg_pol(ll,0,ctheta)
!          Compute for m>0
           onem=one
           do mm=1,ll
             onem=-onem
             work1=ylmcst*sqrt(fact)*onem*ass_leg_pol(ll,mm,ctheta)*sqrt(2._dp)
             ylmr(l0+mm,inpt)=work1*rphase(mm)
             ylmr(l0-mm,inpt)=work1*iphase(mm)
             if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
           end do ! End loop over m
         end do  ! End loop over l
       end if

!      COMPUTE dYlm/dRi
       if (compute_ylmgr) then
!        Loop over angular momentum l
         do ilang=2,mpsang
           ll=ilang-1
           l0=ll**2+ll+1
           fact=1._dp/real(ll*(ll+1),dp)
           ylmcst=sqrt(real(2*ll+1,dp)/four_pi)/rnorm
!          Special case m=0
           work1=ylmcst*plm_dtheta(ll,0,ctheta)
           ylmr_gr(1:3,l0,inpt)=work1*dtheta(1:3)
!          Compute for m>0
           onem=one
           do mm=1,ll
             onem=-onem
             work1=ylmcst*sqrt(fact)*onem*plm_dtheta(ll,mm,ctheta)*sqrt(2._dp)
             work2=ylmcst*sqrt(fact)*onem*plm_dphi  (ll,mm,ctheta)*sqrt(2._dp)
             ylmr_gr(1:3,l0+mm,inpt)=rphase(mm)*work1*dtheta(1:3)-iphase(mm)*work2*dphi(1:3)
             ylmr_gr(1:3,l0-mm,inpt)=iphase(mm)*work1*dtheta(1:3)+rphase(mm)*work2*dphi(1:3)
             if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
           end do ! End loop over m
         end do  ! End loop over l
       end if

!      COMPUTE d2Ylm/dRidRj
       if (compute_ylm2gr) then
         LIBPAW_ALLOCATE(blm,(5,mpsang*mpsang))
         call plm_coeff(blm,mpsang,ctheta)

!        Loop over angular momentum l
         do ilang=2,mpsang
           ll=ilang-1
           l0=ll**2+ll+1
           fact=1._dp/real(ll*(ll+1),dp)
           ylmcst=sqrt(real(2*ll+1,dp)/four_pi)/(rnorm**2)
!          Special case m=0
           ylmr_gr(4,l0,inpt)=ylmcst*(-blm(3,l0)*sphi*sphi+blm(4,l0)*cphi*cphi)
           ylmr_gr(5,l0,inpt)=ylmcst*(-blm(3,l0)*cphi*cphi+blm(4,l0)*sphi*sphi)
           ylmr_gr(6,l0,inpt)=ylmcst*blm(1,l0)
           ylmr_gr(7,l0,inpt)=ylmcst*blm(2,l0)*sphi
           ylmr_gr(8,l0,inpt)=ylmcst*blm(2,l0)*cphi
           ylmr_gr(9,l0,inpt)=ylmcst*(blm(3,l0)+blm(4,l0))*sphi*cphi
!          Compute for m>0
           onem=one
           do mm=1,ll
             onem=-onem;ylmcst2=ylmcst*sqrt(fact)*sqrt(two)
             ylmr_gr(4,l0+mm,inpt)=ylmcst2*((-blm(3,l0+mm)*sphi*sphi+blm(4,l0+mm)*cphi*cphi)*rphase(mm)-&
&             blm(5,l0+mm)*2.d0*cphi*sphi*mm*iphase(mm))
             ylmr_gr(4,l0-mm,inpt)=ylmcst2*((-blm(3,l0+mm)*sphi*sphi+blm(4,l0+mm)*cphi*cphi)*iphase(mm)+&
&             blm(5,l0+mm)*2.d0*cphi*sphi*mm*rphase(mm))
             ylmr_gr(5,l0+mm,inpt)=ylmcst2*((-blm(3,l0+mm)*cphi*cphi+blm(4,l0+mm)*sphi*sphi)*rphase(mm)+&
&             blm(5,l0+mm)*2.d0*cphi*sphi*mm*iphase(mm))
             ylmr_gr(5,l0-mm,inpt)=ylmcst2*((-blm(3,l0+mm)*cphi*cphi+blm(4,l0+mm)*sphi*sphi)*iphase(mm)-&
&             blm(5,l0+mm)*2.d0*cphi*sphi*mm*rphase(mm))
             ylmr_gr(6,l0+mm,inpt)=ylmcst2*blm(1,l0+mm)*rphase(mm)
             ylmr_gr(6,l0-mm,inpt)=ylmcst2*blm(1,l0+mm)*iphase(mm)
             ylmr_gr(7,l0+mm,inpt)=ylmcst2*(blm(2,l0+mm)*sphi*rphase(mm)+&
&             mm*iphase(mm)*cphi*onem*plm_dtheta(ll,mm,ctheta))
             ylmr_gr(7,l0-mm,inpt)=ylmcst2*(blm(2,l0+mm)*sphi*iphase(mm)-&
&             mm*rphase(mm)*cphi*onem*plm_dtheta(ll,mm,ctheta))
             ylmr_gr(8,l0+mm,inpt)=ylmcst2*(blm(2,l0+mm)*cphi*rphase(mm)-&
&             mm*iphase(mm)*sphi*onem*plm_dtheta(ll,mm,ctheta))
             ylmr_gr(8,l0-mm,inpt)=ylmcst2*(blm(2,l0+mm)*cphi*iphase(mm)+&
&             mm*rphase(mm)*sphi*onem*plm_dtheta(ll,mm,ctheta))
             ylmr_gr(9,l0+mm,inpt)=ylmcst2*((blm(3,l0+mm)+blm(4,l0+mm))*sphi*cphi*rphase(mm)-&
&             blm(5,l0+mm)*(sphi*sphi-cphi*cphi)*mm*iphase(mm))
             ylmr_gr(9,l0-mm,inpt)=ylmcst2*((blm(3,l0+mm)+blm(4,l0+mm))*sphi*cphi*iphase(mm)+&
&             blm(5,l0+mm)*(sphi*sphi-cphi*cphi)*mm*rphase(mm))
             if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
           end do ! End loop over m
         end do  ! End loop over l
         LIBPAW_DEALLOCATE(blm)
       end if

!      End condition r<>0
     end if

!    End loop over rr
   end do

!  End condition l<>0
 end if

end subroutine initylmr
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/ys
!! NAME
!! ys
!!
!! FUNCTION
!!  Computes the matrix element <Yl'm'|Slm>
!!
!! INPUTS
!!  integer :: l',m',l,m
!!
!! OUTPUT
!!  complex(dpc) :: ys_val 
!! 
!! NOTES
!! Ylm is the standard complex-valued spherical harmonic, Slm is the real spherical harmonic
!! used througout abinit. <Yl'm'|Slm> is their overlap.
!!
!! PARENTS
!!      m_epjdos,m_paw_sphharm
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine ys(lp,mp,ll,mm,ys_val)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ys'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll,lp,mm,mp
 complex(dpc),intent(out) :: ys_val

!Local variables ---------------------------------------
!scalars
 complex(dpc) :: dmpmm,dmpmmm,m1mm

! *********************************************************************

 
 ys_val = czero

 if(lp==ll .AND. (mp==mm .OR. mp==-mm) ) then
  ! (-1)**mm
   m1mm=cone; if(abs(mod(mm,2))==1) m1mm=-m1mm
  
  ! delta(mp,mm)
   dmpmm=czero; if(mp==mm) dmpmm=cone
  
  ! delta(mp,-mm)
   dmpmmm=czero; if(mp==-mm) dmpmmm=cone

   select case (mm)
       case (0) ! case for S_l0
         ys_val = dmpmm
       case (:-1) ! case for S_lm with m < 0
         ys_val = -(zero,one)*m1mm*sqrthalf*(dmpmmm-m1mm*dmpmm)
       case (1:) ! case for S_lm with m > 0
         ys_val = m1mm*sqrthalf*(dmpmm+m1mm*dmpmmm)
   end select

 end if
 
end subroutine ys
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/lxyz.F90
!! NAME
!! lxyz
!!
!! FUNCTION
!! Computes the matrix element <Yl'm'|L_idir|Ylm>
!!
!! INPUTS
!!   integer :: lp,mp,idir,ll,mm
!!
!! OUTPUT
!!   complex(dpc) :: lidir
!! 
!! NOTES
!!  Ylm is the standard complex-valued spherical harmonic,
!!  idir is the direction in space of L
!!
!! PARENTS
!!      m_paw_sphharm
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine lxyz(lp,mp,idir,ll,mm,lidir)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lxyz'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: idir,ll,lp,mm,mp
 complex(dpc),intent(out) :: lidir

!Local variables ---------------------------------------
!scalars
 complex(dpc) :: jmme, jpme

! *********************************************************************

 jpme=czero; jmme=czero
 if (lp==ll) then
   if (mp==mm+1) jpme=cone*sqrt((ll-mm)*(ll+mm+one))
   if (mp==mm-1) jmme=cone*sqrt((ll-mm+one)*(ll+mm))
 end if
 
 lidir = czero
 if (lp == ll) then 
   select case (idir)
     case (1) ! Lx
       lidir = cone*half*(jpme+jmme)
     case (2) ! Ly
       lidir = -(zero,one)*half*(jpme-jmme)
     case (3) ! Lz
       if (mp == mm) lidir = mm*cone
   end select
 end if

end subroutine lxyz
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/slxyzs
!! NAME
!! slxyzs
!!
!! FUNCTION
!! computes the matrix element <Sl'm'|L_idir|Slm>
!!
!! INPUTS
!!   integer :: lp,mp,idir,ll,mm
!!
!! OUTPUT
!!   complex(dpc) :: sls_val
!! 
!! NOTES
!! Slm is the real spherical harmonic used througout abinit,
!! L_idir is a component of the angular momentum operator.
!! The subroutine computes <S_l'm'|L_idir|S_lm>
!!
!! PARENTS
!!      m_pawdij
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine slxyzs(lp,mp,idir,ll,mm,sls_val)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'slxyzs'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: idir,ll,lp,mm,mp
 complex(dpc),intent(out) :: sls_val

!Local variables ---------------------------------------
!scalars
 integer :: lpp,lppp,mpp,mppp
 complex(dpc) :: lidir,sy_val,ys_val

! *********************************************************************

 sls_val = czero
 
 if (lp == ll) then
   lpp  = ll
   lppp = ll
   do mpp = -lpp, lpp
     call ys(lpp,mpp,lp,mp,sy_val)
     do mppp = -lppp, lppp
       call lxyz(lpp,mpp,idir,lppp,mppp,lidir)
       call ys(lppp,mppp,ll,mm,ys_val)
       sls_val = sls_val + conjg(sy_val)*lidir*ys_val
     end do
   end do
 end if
 
end subroutine slxyzs
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/plm_coeff
!! NAME
!! plm_coeff
!!
!! FUNCTION
!! Compute coefficients depending on Plm and its derivatives where P_lm is a legendre polynomial.
!! They are used to compute the second derivatives of spherical harmonics
!!
!! INPUTS
!!  mpsang=1+ maximum l quantum number
!!  xx= input value
!!
!! OUTPUT
!!  blm(5,mpsang*mpsang)=coefficients depending on Plm and its derivatives where P_lm is a legendre polynome
!!
!! PARENTS
!!      initylmg,m_paw_sphharm
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine plm_coeff(blm,mpsang,xx)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'plm_coeff'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: mpsang
 real(dp),intent(in) :: xx
!arrays
 real(dp),intent(out) :: blm(5,mpsang*mpsang)

!Local variables ---------------------------------------
!scalars
 integer :: il,ilm,ilm0,ilm1,im
 real(dp) :: dplm_dt,d2plm_dt2,llp1,onemx2,plm,sqrx,xsqrx,xx2,yy
 logical :: is_one
 character(len=500) :: msg
!arrays
 real(dp) :: pl_d2(mpsang),plm_d2t(mpsang*mpsang)

!************************************************************************

 if (abs(xx).gt.1.d0) then
   msg = ' plm_coeff :  xx > 1 !'
   MSG_ERROR(msg)
 end if

 blm=zero
 is_one=(abs(abs(xx)-one)<=tol12)
 xx2=xx**2
 onemx2=abs(one-xx2)
 sqrx=sqrt(onemx2)
 xsqrx=xx*sqrt(onemx2)

 call plm_d2theta(mpsang,plm_d2t,xx)
 if (is_one) then
   yy=sign(one,xx)
   call pl_deriv(mpsang,pl_d2,yy)
 end if

 do il=0,mpsang-1
   llp1=dble(il*(il+1))
   ilm0=il*il+il+1
   do im=0,il
     ilm=ilm0+im;ilm1=ilm0-im

     plm      =(-1)**im*ass_leg_pol(il,im,xx)
     dplm_dt  =(-1)**im*plm_dtheta(il,im,xx)
     d2plm_dt2=         plm_d2t(ilm)

     blm(1,ilm)=         two*xsqrx    *dplm_dt+onemx2*d2plm_dt2
     blm(2,ilm)=         (one-two*xx2)*dplm_dt-xsqrx *d2plm_dt2
     blm(3,ilm)=llp1*plm+                             d2plm_dt2
     blm(4,ilm)=        -two*xsqrx    *dplm_dt+xx2   *d2plm_dt2


     if (is_one) then
       if (im==1) then
         blm(5,ilm)=llp1*plm+d2plm_dt2
       end if
       if (im==2) then
         blm(5,ilm)=d2plm_dt2-three*pl_d2(il+1)
       end if
     else
       if(im>0) then
         blm(5,ilm)=plm/onemx2-dplm_dt*xx/sqrx
       end if
     end if

     if (im>0) then
       blm(1,ilm1)=blm(1,ilm)
       blm(2,ilm1)=blm(2,ilm)
       blm(3,ilm1)=blm(3,ilm)
       blm(4,ilm1)=blm(4,ilm)
       blm(5,ilm1)=blm(5,ilm)
     end if

   end do
 end do

end subroutine plm_coeff
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/ass_leg_pol
!! NAME
!! ass_leg_pol
!!
!! FUNCTION
!! Compute the associated Legendre Polynomial Plm(x),
!! using a stable recursion formula.
!! Here m and l are integers satisfying 0<=m<=l,
!! while x lies in the range -1<=x<=1
!!
!! INPUTS
!!  l,m= l,m numbers
!!  xarg=argument of the polynom
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ass_leg_pol(l,m,xarg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ass_leg_pol'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) ::  l,m
 real(dp), intent(in) :: xarg
 real(dp) :: ass_leg_pol

!Local variables-------------------------------
!scalars
 integer :: i,ll
 real(dp) :: pll,polmm,tmp1,sqrx,x
 character(len=100) :: msg

! *************************************************************************

 x=xarg
 if (m.lt.0.or.m.gt.l.or.abs(x).gt.1.d0) then
   if (m.lt.0.or.m.gt.l.or.abs(x).gt.1.d0+1.d-10) then
    msg='Bad choice of l, m or x !'
    MSG_BUG(msg)
   endif
   x=1.d0
 endif

 polmm=1.d0
 if (m>0) then
  sqrx=sqrt(abs((1.d0-x)*(1.d0+x)))
  do i=1,m
   polmm=polmm*(1.0d0-2.0d0*i)*sqrx
  enddo
 endif

 if (l==m) then
  ass_leg_pol=polmm
 else
  tmp1=x*(2.0d0*m+1.0d0)*polmm
  if (l==(m+1)) then
   ass_leg_pol=tmp1
  else
   do ll=m+2,l
    pll=(x*(2.0d0*ll-1.0d0)*tmp1-(ll+m-1.0d0)*polmm)/dble(ll-m)
    polmm=tmp1
    tmp1=pll
   enddo
   ass_leg_pol=pll
  endif
 endif

end function ass_leg_pol
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/plm_d2theta
!! NAME
!! plm_d2theta
!!
!! FUNCTION
!! Compute d2(Plm (cos(theta)))/d(theta)2  where P_lm is a legendre polynome
!!
!! INPUTS
!!  mpsang=1+ maximum l quantum number
!!  xx= input value
!!
!! OUTPUT
!!  plm_d2t(mpsang*mpsang)
!!
!! PARENTS
!!      m_paw_sphharm
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine plm_d2theta(mpsang,plm_d2t,xx)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'plm_d2theta'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: mpsang
 real(dp),intent(in) :: xx
!arrays
 real(dp),intent(out) :: plm_d2t(mpsang*mpsang)

!Local variables ---------------------------------------
!scalars
 integer :: il,ilm,ilmm1,ilmm2,im
 real(dp) :: sqrx
 character(len=500) :: msg

!************************************************************************
 if (abs(xx).gt.1.d0) then
   msg = 'plm_d2theta : xx > 1 !'
   MSG_ERROR(msg)
 end if

 plm_d2t=zero
 if (mpsang>1) then
   sqrx=sqrt(abs((1.d0-xx)*(1.d0+xx)))

   do il=1,mpsang-1
     ilm=il*il+2*il+1
     ilmm1=(il-1)*(il-1)+2*(il-1)+1
!    terme d2(Pll)/dtet2
     plm_d2t(ilm)=(2*il-1)*(sqrx*(plm_d2t(ilmm1)-(-1)**(il-1)*ass_leg_pol(il-1,il-1,xx))+&
&     2.d0*xx*(-1)**(il-1)*plm_dtheta(il-1,il-1,xx))
     plm_d2t(ilm-2*il)=plm_d2t(ilm)
!    terme d2(Pl(l-1))/dtet2
     plm_d2t(ilm-1)=(2*il-1)*(xx*(plm_d2t(ilmm1)-(-1)**(il-1)*ass_leg_pol(il-1,il-1,xx))-&
&     2.d0*sqrx*(-1)**(il-1)*plm_dtheta(il-1,il-1,xx))
     if(il>1) plm_d2t(il*il+2)=plm_d2t(ilm-1)
   end do
!  terme d2(Plm)/dtet2
   if(mpsang>2) then
     do il=2,mpsang-1
       do im=0,il-2
         ilm=il*il+il+1+im
         ilmm1=(il-1)*(il-1)+il+im
         ilmm2=(il-2)*(il-2)+il-1+im
         plm_d2t(ilm)=dble(2*il-1)/dble(il-im)*(xx*(plm_d2t(ilmm1)-(-1)**im*ass_leg_pol(il-1,im,xx))-&
&         2.d0*sqrx*(-1)**im*plm_dtheta(il-1,im,xx))-&
&         dble(il+im-1)/dble(il-im)*plm_d2t(ilmm2)
         plm_d2t(ilm-2*im)=plm_d2t(ilm)
       end do
     end do
   end if
 end if

end subroutine plm_d2theta
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/plm_dphi
!! NAME
!! plm_dphi
!!
!! FUNCTION
!! Compute  m*P_lm(x)/sqrt((1-x^2)where P_lm is a legendre polynome
!!
!! INPUTS
!!  ll= l quantum number
!!  mm= m quantum number
!!  xx= input value
!!
!! OUTPUT
!!  plm_dphi(xx)
!!
!! NOTES
!!  This routine comes from Function Der_Phi_P(L,m,x)
!!  (pwpaw code from N. Holzwarth, implemented by Y. Abraham))
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function plm_dphi(ll,mm,xx)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'plm_dphi'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll,mm
 real(dp) :: plm_dphi
 real(dp),intent(in) :: xx

!Local variables ---------------------------------------
!scalars
 integer :: il,im
 real(dp) :: dosomx2,fact,pll,pmm,pmmp1,somx2
 character(len=500) :: msg

! *********************************************************************

 if (mm.lt.0.or.mm.gt.ll.or.abs(xx).gt.1.d0) then
   msg = 'plm_dphi : mm < 0 or mm > ll or xx > 1 !'
   MSG_ERROR(msg)
 end if

 plm_dphi=zero
 if (mm==0) return

 pmm=one
 dosomx2=one
 if (mm > 0) then
   somx2=sqrt((1-xx)*(1+xx))
   fact=one
   do im=1,mm
     pmm=-pmm*fact
     fact=fact+2
   end do
   if (mm > 1) then
     do im=2,mm
       dosomx2=somx2*dosomx2
     end do
   end if
   pmm=pmm*dosomx2 !due to one more term (-1^M)
 end if
 if(ll==mm) then
   plm_dphi=pmm*mm
 else
   pmmp1=xx*(2*mm+1)*pmm
   if(ll==mm+1) then
     plm_dphi=pmmp1*mm
   else if(ll>=mm+2) then
     do il=mm+2,ll
       pll=(xx*(2*il-1)*pmmp1-(il+mm-1)*pmm)/(il-mm)
       pmm=pmmp1
       pmmp1=pll
     end do
     plm_dphi=pll*mm
   end if
 end if

end function plm_dphi
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/plm_dtheta
!! NAME
!! plm_dtheta
!!
!! FUNCTION
!! Compute -(1-x^2)^1/2*d/dx{P_lm(x)} where P_lm is a legendre polynome
!!
!! INPUTS
!!  ll= l quantum number
!!  mm= m quantum number
!!  xx= input value
!!
!! OUTPUT
!!  plm_dtheta(xx)
!!
!! NOTES
!!  This routine comes from Function Der_Theta_P(L,m,x)
!!  (pwpaw code from N. Holzwarth, implemented by Y. Abraham))
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function plm_dtheta(ll,mm,xx)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'plm_dtheta'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll,mm
 real(dp) :: plm_dtheta
 real(dp),intent(in) :: xx

!Local variables ---------------------------------------
!scalars
 integer :: il,im
 real(dp) :: dosomx2,dpll,dpmm,dpmmp1,fact,pll,pmm,pmmp1,somx2
 character(len=500) :: msg

! *********************************************************************

 if (mm.lt.0.or.mm.gt.ll.or.abs(xx).gt.1.d0) then
   msg = 'plm_dtheta : mm < 0 or mm > ll or xx > 1 !'
   MSG_ERROR(msg)
 end if

 plm_dtheta=zero
 pmm=one
 dpmm=one
 dosomx2=one
 somx2=sqrt((1-xx)*(1+xx))
 if(mm==0)then
   dpmm=zero
 elseif (mm > 0) then
   fact=one
   do im=1,mm
     pmm=-pmm*fact*somx2
     dpmm=-dpmm*fact
     fact=fact+2
   end do
   if(mm>1)then
     do im=2,mm
       dosomx2=dosomx2*somx2
     end do
   end if
   dpmm= dpmm*mm*xx*dosomx2
 end if
 if(ll==mm)then
   plm_dtheta=dpmm
 else
   pmmp1=xx*(2*mm+1)*pmm
   dpmmp1=-(2*mm+1)*somx2*pmm+xx*(2*mm+1)*dpmm
   if(ll==mm+1) then
     plm_dtheta=dpmmp1
   else if(ll>=mm+2)then
     do il=mm+2,ll
       pll=(xx*(2*il-1)*pmmp1-(il+mm-1)*pmm)/(il-mm)
       dpll=(-somx2*(2*il-1)*pmmp1+(xx*(2*il-1)*dpmmp1-(il+mm-1)*dpmm))/(il-mm)
       pmm=pmmp1
       pmmp1=pll
       dpmm=dpmmp1
       dpmmp1=dpll
     end do
     plm_dtheta=dpll
   end if
 end if

end function plm_dtheta
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/pl_deriv
!! NAME
!! pl_deriv
!!
!! FUNCTION
!! Compute d2(Pl (x)))/d(x)2  where P_l is a legendre polynomial
!!
!! INPUTS
!!  mpsang=1+ maximum l quantum number
!!  xx= input value
!!
!! OUTPUT
!!  pl_d2(mpsang*mpsang)
!!
!! PARENTS
!!      m_paw_sphharm
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine pl_deriv(mpsang,pl_d2,xx)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pl_deriv'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: mpsang
 real(dp),intent(in) :: xx
!arrays
 real(dp),intent(out) :: pl_d2(mpsang)

!Local variables ---------------------------------------
!scalars
 integer :: il,ilm
 real(dp) :: il_,il_m1,il_2m1
 character(len=500) :: msg
!arrays
 real(dp) :: pl(mpsang),pl_d1(mpsang)

! *********************************************************************

 if (abs(xx).gt.1.d0) then
   msg = 'pl_deriv : xx > 1 !'
   MSG_ERROR(msg)
 end if

 pl_d2=zero; pl_d1=zero; pl=zero
 pl(1)=one; pl(2)=xx
 pl_d1(1)=zero; pl_d1(2)=one
 pl_d2(1)=zero; pl_d2(2)=zero
 if (mpsang>2) then
   do il=2,mpsang-1
     il_=dble(il);il_m1=dble(il-1);il_2m1=dble(2*il-1)
     ilm=il+1
     pl(ilm)=(il_2m1*xx*pl(ilm-1)-il_m1*pl(ilm-2))/il_
     pl_d1(ilm)=(il_2m1*(xx*pl_d1(ilm-1)+pl(ilm-1))-il_m1*pl_d1(ilm-2))/il_
     pl_d2(ilm)=(il_2m1*(xx*pl_d2(ilm-1)+two*pl_d1(ilm-1))-il_m1*pl_d2(ilm-2))/il_
   end do
 end if

end subroutine pl_deriv
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/mat_mlms2jmj
!! NAME
!! mat_mlms2jmj
!!
!! FUNCTION
!! For a given angular momentum lcor, change a matrix of dimension 2(2*lcor+1)
!! from the Ylm basis to the J,M_J basis if option==1
!!
!! INPUTS
!!  lcor= angular momentum
!!  ndij= ndij = 4
!!  option=  1 matrix in |l,s,m_l,m_s> basis is changed into |l,s,j,m_j> basis
!!           2 matrix in |l,s,j,m_j> basis is changed into |l,s,m_l,m_s> basis
!!  optspin=  1  Spin up are first
!!            2  Spin dn are first
!!  prtvol=printing volume
!!  unitfi=printing file unit ; -1 for no printing
!!  wrt_mode=printing mode in parallel ('COLL' or 'PERS')
!!
!! SIDE EFFECTS
!!  mat_mlms= Input/Ouput matrix in the Ylm basis, size of the matrix is (2*lcor+1,2*lcor+1,ndij)
!!  mat_jmj= Input/Output matrix in the J,M_J basis, size is 2*(2*lcor+1),2*(2*lcor+1)
!!
!! NOTES
!!  usefull only in ndij==4
!!
!! PARENTS
!!      m_pawang,pawprt,setnoccmmp
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine mat_mlms2jmj(lcor,mat_mlms,mat_jmj,ndij,option,optspin,prtvol,unitfi,wrt_mode)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mat_mlms2jmj'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ndij,lcor,option,optspin,prtvol,unitfi
 character(len=4),intent(in) :: wrt_mode
!arrays
 complex(dpc),intent(inout) :: mat_mlms(2*lcor+1,2*lcor+1,ndij)
 complex(dpc),intent(inout) :: mat_jmj(2*(2*lcor+1),2*(2*lcor+1))

!Local variables ---------------------------------------
!scalars
 integer :: ii,im,im1,im2,ispden,jc1,jc2,jj,jm,ll,ml1,ml2,ms1,ms2
 real(dp),parameter :: invsqrt2=one/sqrt2
 real(dp) :: invsqrt2lp1,xj,xmj
 complex(dpc) :: mat_tmp,tmp2
 character(len=9),parameter :: dspinold(6)=(/"up       ","down     ","up-up    ","down-down","up-dn    ","dn-up    "/)
 character(len=9),parameter :: dspin(6)=(/"dn       ","up       ","dn-dn    ","up-up    ","dn-up    ","up-dn    "/)
 character(len=500) :: msg
!arrays
 integer, allocatable :: ind_msml(:,:)
 complex(dpc),allocatable :: mat_mlms2(:,:),mlms2jmj(:,:)

!*********************************************************************

 if(ndij/=4) then
   msg=" ndij/=4 !"
   MSG_BUG(msg)
 end if
 if (option/=1.and.option/=2) then
   msg=' option=/1 and =/2 !'
   MSG_BUG(msg)
 end if
 if (optspin/=1.and.optspin/=2) then
   msg=' optspin=/1 and =/2 !'
   MSG_BUG(msg)
 end if
 
 if (unitfi/=-1) then
   if(option==1) then
     write(msg,'(3a)') ch10,&
&     "matrix in |l,s,m_l,m_s> basis is changed into |l,s,j,m_j> basis"
     call wrtout(unitfi,msg,wrt_mode)
   else if(option==2) then
     write(msg,'(3a)') ch10,&
&     "matrix in |l,s,j,m_j> basis is changed into |l,s,m_l,m_s> basis"
     call wrtout(unitfi,msg,wrt_mode)
   end if
 end if
 
 if(option==1) then
   if(optspin==2) then
     if(abs(prtvol)>2.and.unitfi/=-1)&
&     write(msg,'(3a)') ch10,"assume spin dn is the first in the array"
   else if (optspin==1) then
     if(abs(prtvol)>2.and.unitfi/=-1)&
&     write(msg,'(3a)') ch10,"change array in order that spin dn is the first in the array"
     do ii=1,2*lcor+1
       do jj=1,2*lcor+1
         mat_tmp=mat_mlms(ii,jj,2)
         mat_mlms(ii,jj,2)=mat_mlms(ii,jj,1)
         mat_mlms(ii,jj,1)=mat_tmp
         mat_tmp=mat_mlms(ii,jj,4)
         mat_mlms(ii,jj,4)=mat_mlms(ii,jj,3)
         mat_mlms(ii,jj,3)=mat_tmp
       end do
     end do
!    mat_tmp(:,:,1)=mat_mlms(:,:,2);mat_tmp(:,:,2)=mat_mlms(:,:,1)
!    mat_tmp(:,:,3)=mat_mlms(:,:,4);mat_tmp(:,:,4)=mat_mlms(:,:,3)
!    mat_mlms(:,:,:)=mat_tmp(:,:,:)
   end if
   if(abs(prtvol)>2.and.unitfi/=-1) then
     call wrtout(unitfi,msg,wrt_mode)
   end if
 end if

 if(option==1.and.abs(prtvol)>2.and.unitfi/=-1) then
   do ispden=1,ndij
     write(msg,'(3a)') ch10,&
&     "Input matrix in the Ylm basis for component ",trim(dspin(ispden+2*(ndij/4)))
     call wrtout(unitfi,msg,wrt_mode)
     do im1=1,lcor*2+1
       write(msg,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))')&
&       (mat_mlms(im1,im2,ispden),im2=1,lcor*2+1)
       call wrtout(unitfi,msg,wrt_mode)
     end do
   end do
 end if ! option==1

!--------------- Built indices + allocations
 ll=lcor
 LIBPAW_ALLOCATE(mlms2jmj,(2*(2*ll+1),2*(2*ll+1)))
 mlms2jmj=czero
 LIBPAW_BOUND2_ALLOCATE(ind_msml,BOUNDS(1,2),BOUNDS(-ll,ll))
 LIBPAW_ALLOCATE(mat_mlms2,(2*(2*lcor+1),2*(2*lcor+1)))
 mlms2jmj=czero
 jc1=0
 do ms1=1,2
   do ml1=-ll,ll
     jc1=jc1+1
     ind_msml(ms1,ml1)=jc1
   end do
 end do
!--------------- Change representation of input matrix for ndij==4
 if(option==1) then
   jc1=0
   do ms1=1,2
     do ml1=1,2*ll+1
       jc1=jc1+1
       jc2=0
       do ms2=1,2
         do ml2=1,2*ll+1
           jc2=jc2+1
           if(ms1==ms2) mat_mlms2(jc1,jc2)=mat_mlms(ml1,ml2,ms1)
           if(ms1<ms2) mat_mlms2(jc1,jc2)=mat_mlms(ml1,ml2,3)
           if(ms1>ms2) mat_mlms2(jc1,jc2)=mat_mlms(ml1,ml2,4)
         end do
       end do
     end do
   end do
   if(abs(prtvol)>1.and.unitfi/=-1) then
     write(msg,'(3a)') ch10,"Input matrix in the lms basis for all component"
     call wrtout(unitfi,msg,wrt_mode)
     do im1=1,2*(lcor*2+1)
       write(msg,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))')&
&       (mat_mlms2(im1,im2),im2=1,2*(lcor*2+1))
       call wrtout(unitfi,msg,wrt_mode)
     end do
   end if
 end if  ! option==1

!--------------- built mlms2jmj
!do jj=ll,ll+1    ! the physical value of j are ll-0.5,ll+0.5
!xj(jj)=jj-0.5
 if(ll==0)then
   msg=' ll should not be equal to zero !'
   MSG_BUG(msg)
 end if
 jc1=0
 invsqrt2lp1=one/sqrt(float(2*lcor+1))
 do jj=ll,ll+1
   xj=float(jj)-half
   do jm=-jj,jj-1
     xmj=float(jm)+half
     jc1=jc1+1
     if(nint(xj+0.5)==ll+1) then
       if(nint(xmj+0.5)==ll+1)  then
         mlms2jmj(ind_msml(2,ll),jc1)=1.0   !  J=L+0.5 and m_J=L+0.5
       else if(nint(xmj-0.5)==-ll-1) then
         mlms2jmj(ind_msml(1,-ll),jc1)=1.0   !  J=L+0.5 and m_J=-L-0.5
       else
         mlms2jmj(ind_msml(2,nint(xmj-0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)+xmj+0.5))
         mlms2jmj(ind_msml(1,nint(xmj+0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)-xmj+0.5))
       end if
     end if
     if(nint(xj+0.5)==ll) then
       mlms2jmj(ind_msml(1,nint(xmj+0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)+xmj+0.5))
       mlms2jmj(ind_msml(2,nint(xmj-0.5)),jc1)=-invsqrt2lp1*(sqrt(float(ll)-xmj+0.5))
     end if
   end do
 end do
 if(abs(prtvol)>2.and.unitfi/=-1) then
   write(msg,'(3a)') ch10,"Matrix to go from |M_L,M_S> to |J,M_J>"
   call wrtout(unitfi,msg,wrt_mode)
   do im1=1,2*(lcor*2+1)
     write(msg,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))') (mlms2jmj(im1,im2),im2=1,2*(lcor*2+1))
     call wrtout(unitfi,msg,wrt_mode)
   end do
 end if

 do jm=1,2*(2*ll+1)
   do im=1,2*(2*ll+1)
     tmp2=czero
     do ii=1,2*(2*ll+1)
       do jj=1,2*(2*ll+1)
         if(option==1) then
           tmp2=tmp2+mat_mlms2(ii,jj)*CONJG(mlms2jmj(ii,im))*(mlms2jmj(jj,jm))
         else if(option==2) then
           tmp2=tmp2+mat_jmj(ii,jj)*(mlms2jmj(im,ii))*CONJG(mlms2jmj(jm,jj)) ! inv=t*
         end if
       end do
     end do
     if(option==1) then
       mat_jmj(im,jm)=tmp2
     else if(option==2) then
       mat_mlms2(im,jm)=tmp2
     end if
   end do
 end do
 if(option==1) then
   if (abs(prtvol)>=1.and.unitfi/=-1) then
     write(msg,'(3a)') ch10," Matrix in the J,M_J basis"
     call wrtout(unitfi,msg,wrt_mode)
     do im1=1,2*(lcor*2+1)
       write(msg,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))') (mat_jmj(im1,im2),im2=1,2*(lcor*2+1))
       call wrtout(unitfi,msg,wrt_mode)
     end do
   end if
 else if(option==2) then
   if (abs(prtvol)>=1.and.unitfi/=-1) then
     write(msg,'(3a)') ch10," Matrix in the m_s m_l basis"
     call wrtout(unitfi,msg,wrt_mode)
     do im1=1,2*(lcor*2+1)
       write(msg,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))') (mat_mlms2(im1,im2),im2=1,2*(lcor*2+1))
       call wrtout(unitfi,msg,wrt_mode)
     end do
   end if
   jc1=0
   do ms1=1,2
     do ml1=1,2*ll+1
       jc1=jc1+1
       jc2=0
       do ms2=1,2
         do ml2=1,2*ll+1
           jc2=jc2+1
           if(ms1==ms2) mat_mlms(ml1,ml2,ms1)=mat_mlms2(jc1,jc2)
           if(ms1<ms2) mat_mlms(ml1,ml2,3)=mat_mlms2(jc1,jc2)
           if(ms1>ms2) mat_mlms(ml1,ml2,4)=mat_mlms2(jc1,jc2)
         end do
       end do
     end do
   end do
 end if
 LIBPAW_DEALLOCATE(mlms2jmj)
 LIBPAW_DEALLOCATE(mat_mlms2)
 LIBPAW_DEALLOCATE(ind_msml)

 end subroutine mat_mlms2jmj
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/mat_slm2ylm
!! NAME
!! mat_slm2ylm
!!
!! FUNCTION
!! For a given angular momentum lcor, change a matrix  of dimension (2*lcor+1)
!! from the Slm to the Ylm basis if option==1 or from Ylm to Slm if !option==2
!!
!! INPUTS
!!  lcor= angular momentum, size of the matrix is 2(2*lcor+1)
!!  mat_inp_c= Input matrix
!!  ndij= ndij = 4
!!  option= -1  Change matrix from Slm to Ylm basis
!!           1  Change matrix from Ylm to Slm basis
!!  optspin=  1  Spin up are first
!!            2  Spin dn are first
!!  prtvol=printing volume
!!  unitfi=printing file unit ; -1 for no printing
!!  wrt_mode=printing mode in parallel ('COLL' or 'PERS')
!!
!! OUTPUT
!!  mat_inp_c= Output matrix in Ylm or Slm basis according to option
!!
!! NOTES
!!  usefull only in ndij==4
!!
!! PARENTS
!!      m_pawang,pawprt,setnoccmmp
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine mat_slm2ylm(lcor,mat_inp_c,mat_out_c,ndij,option,optspin,prtvol,unitfi,wrt_mode)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mat_slm2ylm'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ndij,lcor,option,optspin,prtvol,unitfi
 character(len=4),intent(in) :: wrt_mode
!arrays
 complex(dpc) :: mat_inp_c(2*lcor+1,2*lcor+1,ndij),mat_out(2*lcor+1,2*lcor+1,ndij)
 complex(dpc) :: mat_out_c(2*lcor+1,2*lcor+1,ndij)

!Local variables ---------------------------------------
!scalars
 integer :: jm,ii,jj,ll,mm,ispden,im,im1,im2
 real(dp),parameter :: invsqrt2=one/sqrt2
 real(dp) :: onem
 complex(dpc) :: tmp2
 character(len=9),parameter :: dspinc(6)=(/"up       ","down     ","up-up    ","down-down","up-dn    ","dn-up    "/)! optspin 1
 character(len=9),parameter :: dspinc2(6)=(/"up       ","down     ","dn-dn    ","up-up    ","dn-up    ","up-dn    "/)! optspin 2
 character(len=500) :: msg
!arrays
 complex(dpc),allocatable :: slm2ylm(:,:)

! *********************************************************************

 if(ndij/=4) then
   msg=' ndij:=4 !'
   MSG_BUG(msg)
 end if
 if (option/=1.and.option/=2.and.option/=3.and.option/=4) then
   msg=' option=/1 or 2 or 3 or 4 !'
   MSG_BUG(msg)
 end if

 if(abs(prtvol)>2.and.unitfi/=-1) then
   write(msg,'(3a)') ch10, "   mat_slm2ylm"
   call wrtout(unitfi,msg,wrt_mode)
 end if
 
 if(abs(prtvol)>2.and.unitfi/=-1) then
   if(option==1.or.option==3) then
     write(msg,'(3a)') ch10,"matrix in Slm basis is changed into Ylm basis"
     call wrtout(unitfi,msg,wrt_mode)
   else if(option==2.or.option==4) then
     write(msg,'(3a)') ch10,"matrix in Ylm basis is changed into Slm basis"
     call wrtout(unitfi,msg,wrt_mode)
   end if
 end if

 ll=lcor
 LIBPAW_ALLOCATE(slm2ylm,(2*ll+1,2*ll+1))
 slm2ylm=czero
 mat_out=zero
 mat_out_c=czero
 do im=1,2*ll+1
   mm=im-ll-1;jm=-mm+ll+1
   onem=dble((-1)**mm)
   if (mm> 0) then
     slm2ylm(im,im)= cmplx(onem*invsqrt2,zero,kind=dp)
     slm2ylm(jm,im)= cmplx(invsqrt2,     zero,kind=dp)
   end if
   if (mm==0) then
     slm2ylm(im,im)=cone
   end if
   if (mm< 0) then
     slm2ylm(im,im)= cmplx(zero,     invsqrt2,kind=dp)
     slm2ylm(jm,im)=-cmplx(zero,onem*invsqrt2,kind=dp)
   end if
 end do
 if(abs(prtvol)>2.and.unitfi/=-1) then
   do ispden=1,ndij
     if(optspin==1) then
       if(option==1.or.option==3)&
&       write(msg,'(3a)') ch10,&
&       "Input matrix in the Slm basis for component ",trim(dspinc(ispden+2*(ndij/4)))
       if(option==2.or.option==3)&
&       write(msg,'(3a)') ch10,&
&       "Input matrix in the Ylm basis for component ",trim(dspinc(ispden+2*(ndij/4)))
     else
       if(option==1.or.option==3)&
&       write(msg,'(3a)') ch10,&
&       "Input matrix in the Slm basis for component ",trim(dspinc2(ispden+2*(ndij/4)))
       if(option==2.or.option==3)&
&       write(msg,'(3a)') ch10,&
&       "Input matrix in the Ylm basis for component ",trim(dspinc2(ispden+2*(ndij/4)))
     end if
     call wrtout(unitfi,msg,wrt_mode)
     do im1=1,lcor*2+1
       write(msg,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
&       (mat_inp_c(im1,im2,ispden),im2=1,lcor*2+1)
       call wrtout(unitfi,msg,wrt_mode)
     end do
   end do
 end if
 do ispden=1,ndij
   do jm=1,2*ll+1
     do im=1,2*ll+1
       tmp2=czero
       do ii=1,2*ll+1
         do jj=1,2*ll+1
           if(option==1) then
             tmp2=tmp2+mat_inp_c(ii,jj,ispden)*(slm2ylm(im,ii))*CONJG(slm2ylm(jm,jj))
           else if(option==2) then
             tmp2=tmp2+mat_inp_c(ii,jj,ispden)*CONJG(slm2ylm(ii,im))*(slm2ylm(jj,jm))
           end if
         end do
       end do
       mat_out_c(im,jm,ispden)=tmp2
     end do
   end do
 end do ! ispden
 do ii=1,2*ll+1
   do jj=1,2*ll+1
     mat_out(ii,jj,1)=real(mat_out_c(ii,jj,1))
     mat_out(ii,jj,2)=real(mat_out_c(ii,jj,2))
     mat_out(ii,jj,3)=real(mat_out_c(ii,jj,3))
     mat_out(ii,jj,4)=aimag(mat_out_c(ii,jj,3))
!    check that n_{m,m'}^{alpha,beta}=conjg(n_{m',m"}^{beta,alpha}).
     if((abs(aimag(mat_out_c(ii,jj,3))+aimag(mat_out_c(jj,ii,4))).ge.0.0001).or. &
&     (abs(real(mat_out_c(ii,jj,3))-real(mat_out_c(jj,ii,4))).ge.0.0001)) then
       write(msg,'(a,4f10.4)') &
&       ' prb with mat_out_c ',mat_out_c(ii,jj,3),mat_out_c(ii,jj,4)
       MSG_BUG(msg)
     end if
   end do
 end do

 LIBPAW_DEALLOCATE(slm2ylm)

end subroutine mat_slm2ylm
!!***

!----------------------------------------------------------------------

END MODULE m_paw_sphharm
!!***



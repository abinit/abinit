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
!! Copyright (C) 2013-2020 ABINIT group (MT, FJ, NH, TRangel)
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
 public :: ylmc            ! Complex Spherical harmonics for l<=3.
 public :: ylmcd           ! First derivative of complex Ylm wrt theta and phi up to l<=3
 public :: ylm_cmplx       ! All (complex) spherical harmonics for lx<=4
 public :: initylmr        ! Real Spherical Harmonics on a set of vectors
 public :: ys              ! Matrix element <Yl'm'|Slm>
 public :: lxyz            ! Matrix element <Yl'm'|L_idir|Ylm>
 public :: slxyzs          ! Matrix element <Sl'm'|L_idir|Slm>
 public :: plm_coeff       ! Coefficients depending on Plm used to compute the 2nd der of Ylm
 public :: ass_leg_pol     ! Associated Legendre Polynomial Plm(x)
 public :: plm_d2theta     ! d2(Plm (cos(theta)))/d(theta)2 (P_lm= ass. legendre polynomial)
 public :: plm_dphi        ! m*P_lm(x)/sqrt((1-x^2)  (P_lm= ass. legendre polynomial)
 public :: plm_dtheta      ! -(1-x^2)^1/2*d/dx{P_lm(x)} (P_lm= ass. legendre polynomial)
 public :: pl_deriv        ! d2(Pl (x)))/d(x)2  where P_l is a legendre polynomial
 public :: mkeuler         ! For a given symmetry operation, determines the corresponding Euler angles
 public :: dble_factorial  ! Compute factorial of an integer; returns a double precision real
 public :: dbeta           ! Calculate the rotation matrix d^l_{m{\prim}m}(beta)
 public :: phim            ! Computes Phi_m[theta]=Sqrt[2] cos[m theta],      if m>0
                           !                       Sqrt[2] sin[Abs(m) theta], if m<0
                           !                       1                        , if m=0
 public :: mat_mlms2jmj    ! Change a matrix from the Ylm basis to the J,M_J basis
 public :: mat_slm2ylm     ! Change a matrix from the Slm to the Ylm basis or from Ylm to Slm
 public :: create_slm2ylm  ! For a given angular momentum lcor, compute slm2ylm
 public :: create_mlms2jmj ! For a given angular momentum lcor, give the rotation matrix msml2jmj
 public :: setsym_ylm      ! Compute rotation matrices expressed in the basis of real spherical harmonics
 public :: gaunt           ! Gaunt coeffients for complex Yml
 public :: realgaunt       ! Compute "real Gaunt coefficients" with "real spherical harmonics"
 public :: setnabla_ylm    ! Compute rotation matrices expressed in the basis of real spherical harmonics
 public :: nablarealgaunt  ! Compute the integrals of nablaSlimi.nablaySjmj Slkmk on the PAW spheres
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
  new_ylmc = SQRT((2*il+1)*dble_factorial(il-ABS(im))/(dble_factorial(il+ABS(im))*four_pi)) * &
&   ass_leg_pol(il,ABS(im),costh) * CMPLX(cosphi,sinphi)**ABS(im)
  if (im<0) new_ylmc=(-one)**(im)*CONJG(new_ylmc)

  if (ABS(new_ylmc-ylmc)>tol6) then
    !MSG_WARNING("Check new_ylmc")
    !write(std_out,*)"il,im,new_ylmc, ylmc",il,im,new_ylmc,ylmc
    !write(std_out,*)"fact",SQRT((2*il+1)*dble_factorial(il-ABS(im))/(dble_factorial(il+ABS(im))*four_pi))
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
!!
!! SOURCE

subroutine ylmcd(il,im,kcart,dth,dphi)

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
 real(dp) :: cosphi,costh,costhreephi,costwophi,r,rxy,sinphi,sinth,sinthreephi,sintwophi,c
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
     c = SQRT(21.d0/(64.d0*pi))
     dth= -c*      (15.d0*costh**3-11.d0*costh)*            CMPLX(cosphi,sinphi)
     dphi=-c*sinth*( 5.d0*costh**2-1          )*(0.d0,1.d0)*CMPLX(cosphi,sinphi)
   else if (ABS(im)==2) then
     c = SQRT(105.d0/(32.d0*pi))
     dth =c*(2.d0*sinth*costh**2-sinth**3)   *CMPLX(costwophi,sintwophi)
     dphi=c*(2.d0*sinth**2*costh)*(0.d0,1.d0)*CMPLX(costwophi,sintwophi)
   else if (abs(im)==3) then
     dth =-SQRT(35.d0/(64.d0*pi))*3.d0*sinth**2*costh*CMPLX(costhreephi,sinthreephi)
     dphi=-SQRT(35.d0/(64.d0*pi))*sinth**3*(0.d0,3.d0)*CMPLX(costhreephi,sinthreephi)
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
!!      m_mlwfovlp
!!
!! CHILDREN
!!
!! SOURCE

subroutine ylm_cmplx(lx,ylm,xx,yy,zz)

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
!!      m_mlwfovlp,m_paw_finegrid,m_paw_mkrho,m_paw_overlap,m_paw_pwaves_lmn
!!      m_pawang,m_positron,m_rec
!!
!! CHILDREN
!!
!! SOURCE

subroutine initylmr(mpsang,normchoice,npts,nrm,option,rr,ylmr,ylmr_gr)

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
!!
!! SOURCE

subroutine ys(lp,mp,ll,mm,ys_val)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll,lp,mm,mp
 complex(dpc),intent(out) :: ys_val

!Local variables ---------------------------------------
 !scalars
 real(dp) :: d_lp_ll,d_mp_mm,d_mp_mbar,d_mp_am,d_mp_ambar,powm,powam

! *********************************************************************


 ys_val = czero

 d_lp_ll = zero
 if (lp .EQ. ll) d_lp_ll = one

 d_mp_mm = zero
 if (mp .EQ. mm) d_mp_mm = one

 d_mp_mbar = zero
 if (mp .EQ. -mm) d_mp_mbar = one

 d_mp_am = zero
 if (mp .EQ. abs(mm)) d_mp_am = one

 d_mp_ambar = zero
 if (mp .EQ. -abs(mm)) d_mp_ambar = one

 powm=-one
 if (mod(mm,2) .EQ. 0) powm = one

 powam=-one
 if (mod(abs(mm),2) .EQ. 0) powam = one

 select case (mm)
 case (0) ! case for S_l0
    ys_val = cone*d_lp_ll*d_mp_mm
 case (:-1) ! case for S_lm with m < 0
    ys_val = (zero,one)*sqrthalf*powm*d_lp_ll*(-d_mp_am+powam*d_mp_ambar)
 case (1:) ! case for S_lm with m > 0
    ys_val = cone*sqrthalf*d_lp_ll*(powm*d_mp_mm+d_mp_mbar)
 end select

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
!!
!! SOURCE

subroutine lxyz(lp,mp,idir,ll,mm,lidir)

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
!!      m_orbmag,m_paw_nmr,m_pawdij
!!
!! CHILDREN
!!
!! SOURCE

subroutine slxyzs(lp,mp,idir,ll,mm,sls_val)

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
!!      m_initylmg,m_paw_sphharm
!!
!! CHILDREN
!!
!! SOURCE

subroutine plm_coeff(blm,mpsang,xx)

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
!!
!! SOURCE

subroutine plm_d2theta(mpsang,plm_d2t,xx)

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
!!
!! SOURCE

subroutine pl_deriv(mpsang,pl_d2,xx)

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

!!****f* m_paw_sphharm/mkeuler
!! NAME
!! mkeuler
!!
!! FUNCTION
!! For a given symmetry operation, determines the corresponding Euler angles
!!
!! INPUTS
!!  rot(3,3)= symmetry matrix
!!
!! OUTPUT
!!  cosalp=  cos(alpha) with alpha=Euler angle 1
!!  cosbeta= cos(beta)  with beta =Euler angle 2
!!  cosgam=  cos(gamma) with gamma=Euler angle 3
!!  isn= error code (0 if the routine exit normally)
!!  sinalp= sin(alpha) with alpha=Euler angle 1
!!  sinbeta= sin(beta)  with beta=Euler angle 2
!!  singam= sin(gamma) with gamma=Euler angle 3
!!
!! NOTES
!!  This file comes from the file crystal_symmetry.f
!!  by N.A.W. Holzwarth and A. Tackett for the code pwpaw
!!  XG20200718 However, this routine was not accurate in the determination
!!  of beta when cosbeta was close to one (indeed this is a special case). 
!!  This has been corrected. Moreover, sinbeta has been made an output in order
!!  to allow accurate calculations in dbeta. Also, tolerances have been made consistent.
!!
!! PARENTS
!!      m_paw_sphharm
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkeuler(rot,cosbeta,sinbeta,cosalp,sinalp,cosgam,singam,isn)

!Arguments ---------------------------------------------
!scalars
 integer,intent(out) :: isn
 real(dp),intent(out) :: cosalp,cosbeta,cosgam,sinalp,sinbeta,singam
!arrays
 real(dp),intent(in) :: rot(3,3)

!Local variables ---------------------------------------
!scalars
 integer :: ier
 real(dp) :: check,sinbeta2
 character(len=500) :: msg

! *********************************************************************

 do isn= -1,1,2

!Old coding, inaccurate
!  cosbeta=real(isn)*rot(3,3)
!  if(abs(1._dp-cosbeta*cosbeta)<tol10) then
!    sinbeta=zero
!  else
!    sinbeta=sqrt(1._dp-cosbeta*cosbeta)
!  end if
!  if (abs(sinbeta).gt.tol10)  then
!    cosalp=isn*rot(3,1)/sinbeta
!    sinalp=isn*rot(3,2)/sinbeta
!    cosgam=-isn*rot(1,3)/sinbeta
!    singam=isn*rot(2,3)/sinbeta
!  else
!    cosalp=isn*rot(1,1)/cosbeta
!    sinalp=isn*rot(1,2)/cosbeta
!    cosgam=one
!    singam=zero
!  end if

!New coding, more accurate
   cosbeta=real(isn)*rot(3,3)
   sinbeta2=rot(1,3)**2+rot(2,3)**2
   if(sinbeta2<tol8**2)then
     sinbeta=zero
     cosalp=isn*rot(1,1)/cosbeta
     sinalp=isn*rot(1,2)/cosbeta
     cosgam=one
     singam=zero
   else
     sinbeta=sqrt(sinbeta2)
     cosalp=isn*rot(3,1)/sinbeta
     sinalp=isn*rot(3,2)/sinbeta
     cosgam=-isn*rot(1,3)/sinbeta
     singam=isn*rot(2,3)/sinbeta
   end if
!

!  Check matrix:
   ier=0
   check=cosalp*cosbeta*cosgam-sinalp*singam
   if (abs(check-isn*rot(1,1))>tol8) ier=ier+1
   check=sinalp*cosbeta*cosgam+cosalp*singam
   if (abs(check-isn*rot(1,2))>tol8) ier=ier+1
   check=-sinbeta*cosgam
   if (abs(check-isn*rot(1,3))>tol8) ier=ier+1
   check=-cosalp*cosbeta*singam-sinalp*cosgam
   if (abs(check-isn*rot(2,1))>tol8) ier=ier+1
   check=-sinalp*cosbeta*singam+cosalp*cosgam
   if (abs(check-isn*rot(2,2))>tol8) ier=ier+1
   check=sinbeta*singam
   if (abs(check-isn*rot(2,3))>tol8) ier=ier+1
   check=cosalp*sinbeta
   if (abs(check-isn*rot(3,1))>tol8) ier=ier+1
   check=sinalp*sinbeta
   if (abs(check-isn*rot(3,2))>tol8) ier=ier+1
   if (ier.eq.0) return
 end do

 isn=0
 write(msg, '(7a)' )&
& 'Error during determination of symetries!',ch10,&
& 'Action: check your input file:',ch10,&
& 'unit cell vectors and/or atoms positions',ch10,&
& 'have to be given with a better precision.'
 MSG_ERROR(msg)

end subroutine mkeuler
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/dble_factorial
!! NAME
!! dble_factorial
!!
!! FUNCTION
!! PRIVATE function
!! Calculates N! as a double precision real.
!!
!! INPUTS
!!   nn=input integer
!!
!! OUTPUT
!!   factorial= N! (double precision)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental function dble_factorial(nn)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: nn
 real(dp) :: dble_factorial

!Local variables ---------------------------------------
!scalars
 integer :: ii

! *********************************************************************

 dble_factorial=one
 do ii=2,nn
   dble_factorial=dble_factorial*ii
 end do

end function dble_factorial
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/dbeta
!! NAME
!! dbeta
!!
!! FUNCTION
!!  Calculate the rotation matrix d^l_{m{\prim}m}(beta) using Eq. 4.14 of
!!  M.E. Rose, Elementary Theory of Angular Momentum,
!!             John Wiley & Sons, New-York, 1957
!!
!! INPUTS
!!  cosbeta= cosinus of beta (=Euler angle)
!!  sinbeta= sinus of beta (=Euler angle)
!!  ll= index l
!!  mm= index m
!!  mp= index m_prime
!!
!! OUTPUT
!!  dbeta= rotation matrix
!!
!! NOTES
!!  - This file comes from the file crystal_symmetry.f
!!    by N.A.W. Holzwarth and A. Tackett for the code pwpaw
!!  - Assume l relatively small so that factorials do not cause
!!    roundoff error
!!  - XG20200718 This routine was inaccurate when cosbeta was close to one or minus one. 
!!    This has been fixed by adding sinbeta argument obtained from mkeuler. 
!!    Tolerances have been adjusted as well.
!!
!! PARENTS
!!     m_paw_sphharm
!!
!! CHILDREN
!!
!! SOURCE

function dbeta(cosbeta,sinbeta,ll,mp,mm)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll,mm,mp
 real(dp) :: dbeta
 real(dp),intent(in) :: cosbeta,sinbeta

!Local variables ------------------------------
!scalars
 integer,parameter :: mxterms=200
 integer :: ii,ina,inb,inc,ml,ms
 real(dp) :: arg,cosbetab2,pref,sinbetab2,sum,tt

!************************************************************************
 dbeta=zero

!Special cases
 if (abs(cosbeta-1._dp).lt.tol10) then
   if (mp.eq.mm) dbeta=1
 else if (abs(cosbeta+1._dp).lt.tol10) then
   if (mp.eq.-mm) dbeta=(-1)**(ll+mm)
 else
!  General case

!!!!! Old coding
!!  This is inaccurate when cosbeta is close to -1
!   cosbetab2=sqrt((1+cosbeta)*0.5_dp)
!!  This is inaccurate when cosbeta is close to +1
!   sinbetab2=sqrt((1-cosbeta)*0.5_dp)
!!!!! End old coding, begin new coding
  if(cosbeta>-tol8)then
    !If cosbeta is positive, cosbeta2 is positive with value >0.7, so one can divide by cosbetab2
    cosbetab2=sqrt((1+cosbeta)*half)
    sinbetab2=sinbeta*half/cosbetab2
  else
    !If cosbeta is negative, sinbeta2 is positive with value >0.7, so one can divide by sinbetab2
    sinbetab2=sqrt((1-cosbeta)*half)
    cosbetab2=sinbeta*half/sinbetab2
  endif
!!!!! End of new coding

   ml=max(mp,mm)
   ms=min(mp,mm)
   if (ml.ne.mp) sinbetab2=-sinbetab2
   tt=-(sinbetab2/cosbetab2)**2
   pref=sqrt((dble_factorial(ll-ms)*dble_factorial(ll+ml))&
&   /(dble_factorial(ll+ms)*dble_factorial(ll-ml)))&
&   /dble_factorial(ml-ms)*(cosbetab2**(2*ll+ms-ml))&
&   *((-sinbetab2)**(ml-ms))
   sum=1._dp
   arg=1._dp
   ina=ml-ll
   inb=-ms-ll
   inc=ml-ms+1
   do ii=1,mxterms
     if (ina.eq.0.or.inb.eq.0) exit
     arg=(arg*ina*inb*tt)/(ii*inc)
     sum=sum+arg
     ina=ina+1
     inb=inb+1
     inc=inc+1
   end do
   dbeta=pref*sum
 end if

end function dbeta
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/phim
!! NAME
!! phim
!!
!! FUNCTION
!! Computes Phi_m[theta]=Sqrt[2] cos[m theta],      if m>0
!!                       Sqrt[2] sin[Abs(m) theta], if m<0
!!                       1                        , if m=0
!!
!! INPUTS
!!  costeta= cos(theta)  (theta= input angle)
!!  mm = index m
!!  sinteta= sin(theta)  (theta= input angle)
!!
!! OUTPUT
!!  phim= Phi_m(theta) (see above)
!!
!! NOTES
!!  - This file comes from the file crystal_symmetry.f
!!    by N.A.W. Holzwarth and A. Tackett for the code pwpaw
!!
!! PARENTS
!!     m_paw_sphharm
!!
!! CHILDREN
!!
!! SOURCE

pure function phim(costheta,sintheta,mm)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: mm
 real(dp) :: phim
 real(dp),intent(in) :: costheta,sintheta

! *********************************************************************

 if (mm==0)  phim=one
 if (mm==1)  phim=sqrt2*costheta
 if (mm==-1) phim=sqrt2*sintheta
 if (mm==2)  phim=sqrt2*(costheta*costheta-sintheta*sintheta)
 if (mm==-2) phim=sqrt2*two*sintheta*costheta
 if (mm==3)  phim=sqrt2*&
& (costheta*(costheta*costheta-sintheta*sintheta)&
& -sintheta*two*sintheta*costheta)
 if (mm==-3) phim=sqrt2*&
& (sintheta*(costheta*costheta-sintheta*sintheta)&
& +costheta*two*sintheta*costheta)

 end function phim
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
!!  mat_mlms= Input/Output matrix in the Ylm basis, size of the matrix is (2*lcor+1,2*lcor+1,ndij)
!!  mat_jmj= Input/Output matrix in the J,M_J basis, size is 2*(2*lcor+1),2*(2*lcor+1)
!!
!! NOTES
!!  usefull only in ndij==4
!!
!! PARENTS
!!      m_paw_correlations,m_paw_tools,m_pawang
!!
!! CHILDREN
!!
!! SOURCE

subroutine mat_mlms2jmj(lcor,mat_mlms,mat_jmj,ndij,option,optspin,prtvol,unitfi,wrt_mode)

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
!!      m_paw_correlations,m_paw_tools,m_pawang
!!
!! CHILDREN
!!
!! SOURCE

subroutine mat_slm2ylm(lcor,mat_inp_c,mat_out_c,ndij,option,optspin,prtvol,unitfi,wrt_mode)

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

!!****f* m_paw_sphharm/create_slm2ylm
!! NAME
!! create_slm2ylm
!!
!! FUNCTION
!! For a given angular momentum lcor, compute slm2ylm.
!!
!! INPUTS
!!  lcor= angular momentum, size of the matrix is 2(2*lcor+1)
!!
!! OUTPUT
!!  slm2ylm(2lcor+1,2lcor+1) = rotation matrix.
!!
!! NOTES
!!  useful only in ndij==4
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine create_slm2ylm(lcor,slmtwoylm)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lcor
!arrays
 complex(dpc),intent(out) :: slmtwoylm(2*lcor+1,2*lcor+1)

!Local variables ---------------------------------------
!scalars
 integer :: jm,ll,mm,im
 real(dp),parameter :: invsqrt2=one/sqrt2
 real(dp) :: onem
!arrays

! *********************************************************************

 ll=lcor
 slmtwoylm=czero
 do im=1,2*ll+1
   mm=im-ll-1;jm=-mm+ll+1
   onem=dble((-1)**mm)
   if (mm> 0) then
     slmtwoylm(im,im)= cmplx(onem*invsqrt2,zero,kind=dp)
     slmtwoylm(jm,im)= cmplx(invsqrt2,     zero,kind=dp)
   end if
   if (mm==0) then
     slmtwoylm(im,im)=cone
   end if
   if (mm< 0) then
     slmtwoylm(im,im)= cmplx(zero,     invsqrt2,kind=dp)
     slmtwoylm(jm,im)=-cmplx(zero,onem*invsqrt2,kind=dp)
   end if
 end do

end subroutine create_slm2ylm
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/create_mlms2jmj
!! NAME
!! create_mlms2jmj
!!
!! FUNCTION
!! For a given angular momentum lcor, give the rotation matrix msml2jmj
!!
!! INPUTS
!!  lcor= angular momentum
!!
!! SIDE EFFECTS
!!  mlms2jmj= rotation matrix
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine create_mlms2jmj(lcor,mlmstwojmj)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lcor
!arrays
 complex(dpc),intent(out) :: mlmstwojmj(2*(2*lcor+1),2*(2*lcor+1))

!Local variables ---------------------------------------
!scalars
 integer :: jc1,jj,jm,ll,ml1,ms1
 real(dp) :: invsqrt2lp1,xj,xmj
 character(len=500) :: msg
!arrays
 integer, allocatable :: ind_msml(:,:)
 complex(dpc),allocatable :: mat_mlms2(:,:)

!*********************************************************************

!--------------- Built indices + allocations
 ll=lcor
 mlmstwojmj=czero
 LIBPAW_BOUND2_ALLOCATE(ind_msml,BOUNDS(1,2),BOUNDS(-ll,ll))
 LIBPAW_ALLOCATE(mat_mlms2,(2*(2*lcor+1),2*(2*lcor+1)))
 mlmstwojmj=czero
 jc1=0
 do ms1=1,2
   do ml1=-ll,ll
     jc1=jc1+1
     ind_msml(ms1,ml1)=jc1
   end do
 end do

!--------------- built mlmstwojmj
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
         mlmstwojmj(ind_msml(2,ll),jc1)=1.0   !  J=L+0.5 and m_J=L+0.5
       else if(nint(xmj-0.5)==-ll-1) then
         mlmstwojmj(ind_msml(1,-ll),jc1)=1.0   !  J=L+0.5 and m_J=-L-0.5
       else
         mlmstwojmj(ind_msml(2,nint(xmj-0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)+xmj+0.5))
         mlmstwojmj(ind_msml(1,nint(xmj+0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)-xmj+0.5))
       end if
     end if
     if(nint(xj+0.5)==ll) then
       mlmstwojmj(ind_msml(1,nint(xmj+0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)+xmj+0.5))
       mlmstwojmj(ind_msml(2,nint(xmj-0.5)),jc1)=-invsqrt2lp1*(sqrt(float(ll)-xmj+0.5))
     end if
   end do
 end do

 LIBPAW_DEALLOCATE(ind_msml)
 LIBPAW_DEALLOCATE(mat_mlms2)

end subroutine create_mlms2jmj
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/setsym_ylm
!! NAME
!! setsym_ylm
!!
!! FUNCTION
!! Compute rotation matrices expressed in the basis of real spherical harmonics
!! This coefficients are used later to symmetrize PAW on-site quantities (rhoij, dij, ...).
!!
!! INPUTS
!!  gprimd(3,3)==dimensional primitive translations for reciprocal space (bohr^-1)
!!  lmax=value of lmax mentioned at the second line of the psp file
!!  nsym=number of symmetry elements in space group
!!  pawprtvol=control print volume and debugging output
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  sym(3,3,nsym)=symmetries of group in terms of operations on primitive translations
!!
!! OUTPUT
!!  zarot(2*lmax+1,2*lmax+1,lmax+1,nsym)=coefficients of the
!!      transformation of real spherical harmonics
!!      under the symmetry operations
!!
!! NOTES
!!  Typical use: sym(:,:,:) is symrec(:,:,:) (rotations in reciprocal space)
!!               because we need symrel^-1 (=transpose[symrec])
!!               to symmetrize quantities.
!!
!!  - This file comes from the file crystal_symmetry.f
!!    by N.A.W. Holzwarth and A. Tackett for the code pwpaw
!!  - Uses sign & phase convension of  M. E. Rose, Elementary Theory of Angular
!!    Momentum, John Wiley & Sons,. inc. 1957)
!!    zalpha = exp(-i*alpha)   zgamma = exp (-i*gamma)
!!  - Assumes each transformation  can be expressed in terms of 3 Euler
!!    angles with or without inversion
!!
!!  Reference for evaluation of rotation matrices in the basis of real SH:
!!  Blanco M.A., Florez M. and Bermejo M.
!!  Journal of Molecular Structure: THEOCHEM, Volume 419, Number 1, 8 December 1997 , pp. 19-27(9)
!!  http://www.unioviedo.es/qcg/art/Theochem419-19-ov-BF97-rotation-matrices.pdf
!!
!! PARENTS
!!      m_berryphase_new,m_bethe_salpeter,m_dfpt_looppert,m_gstate,m_nonlinear
!!      m_orbmag,m_respfn_driver,m_screening_driver,m_sigma_driver
!!      m_wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine setsym_ylm(gprimd,lmax,nsym,pawprtvol,rprimd,sym,zarot)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lmax,nsym,pawprtvol
!arrays
 integer,intent(in) :: sym(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3)
 real(dp),intent(out) :: zarot(2*lmax+1,2*lmax+1,lmax+1,nsym)

!Local variables ------------------------------
!scalars
 integer :: i1,ii,il,irot,isn,j1,jj,k1,ll,mm,mp
 real(dp) :: cosalp,cosbeta,cosgam,sinalp,sinbeta,singam
 character(len=1000) :: msg
!arrays
 real(dp) :: prod(3,3),rot(3,3)
!************************************************************************

 if (abs(pawprtvol)>=3) then
   write(msg,'(8a,i4)') ch10,&
&   ' PAW TEST:',ch10,&
&   ' ==== setsym_ylm: rotation matrices in the basis ============',ch10,&
&   ' ====              of real spherical harmonics    ============',ch10,&
&   '  > Number of symmetries (nsym)=',nsym
   call wrtout(std_out,msg,'COLL')
 end if

 zarot=zero

 do irot=1,nsym

   if (abs(pawprtvol)>=3) then
     write(msg,'(a,i2,a,9i2,a)') '   >For symmetry ',irot,' (',sym(:,:,irot),')'
     call wrtout(std_out,msg,'COLL')
   end if

!  === l=0 case ===
   zarot(1,1,1,irot)=one

!  === l>0 case ===
   if (lmax>0) then
!    Calculate the rotations in the cartesian basis
     rot=zero;prod=zero
     do k1=1,3
       do j1=1,3
         do i1=1,3
           prod(i1,j1)=prod(i1,j1)+sym(i1,k1,irot)*rprimd(j1,k1)
         end do
       end do
     end do
     do j1=1,3
       do i1=1,3
         do k1=1,3
           rot(i1,j1)=rot(i1,j1)+gprimd(i1,k1)*prod(k1,j1)
         end do
         if(abs(rot(i1,j1))<tol10) rot(i1,j1)=zero
       end do
     end do
     call mkeuler(rot,cosbeta,sinbeta,cosalp,sinalp,cosgam,singam,isn)
     do ll=1,lmax
       il=(isn)**ll
       do mp=-ll,ll
         jj=mp+ll+1
         do mm=-ll,ll
           ii=mm+ll+1

!          Formula (47) from the paper of Blanco et al
           zarot(ii,jj,ll+1,irot)=il&
&           *(phim(cosalp,sinalp,mm)*phim(cosgam,singam,mp)*sign(1,mp)&
           *(dbeta(cosbeta,sinbeta,ll,abs(mp),abs(mm))&
&           +(-1._dp)**mm*dbeta(cosbeta,sinbeta,ll,abs(mm),-abs(mp)))*half&
&           -phim(cosalp,sinalp,-mm)*phim(cosgam,singam,-mp)*sign(1,mm)&
           *(dbeta(cosbeta,sinbeta,ll,abs(mp),abs(mm))&
&           -(-1._dp)**mm*dbeta(cosbeta,sinbeta,ll,abs(mm),-abs(mp)))*half)
         end do
       end do
     end do
   end if   ! lmax case

   if (abs(pawprtvol)>=3) then
     if(lmax>0) then
       write(msg,'(2a,3(3(2x,f7.3),a))') &
&       '    Rotation matrice for l=1:',ch10,&
&       (zarot(1,jj,2,irot),jj=1,3),ch10,&
&       (zarot(2,jj,2,irot),jj=1,3),ch10,&
&       (zarot(3,jj,2,irot),jj=1,3)
       call wrtout(std_out,msg,'COLL')
     end if
     if(lmax>1) then
       write(msg,'(2a,5(5(2x,f7.3),a))') &
&       '    Rotation matrice for l=2:',ch10,&
&       (zarot(1,jj,3,irot),jj=1,5),ch10,&
&       (zarot(2,jj,3,irot),jj=1,5),ch10,&
&       (zarot(3,jj,3,irot),jj=1,5),ch10,&
&       (zarot(4,jj,3,irot),jj=1,5),ch10,&
&       (zarot(5,jj,3,irot),jj=1,5)
       call wrtout(std_out,msg,'COLL')
     end if
     if(lmax>2) then
       write(msg,'(2a,7(7(2x,f7.3),a))') &
&       '    Rotation matrice for l=3:',ch10,&
&       (zarot(1,jj,4,irot),jj=1,7),ch10,&
&       (zarot(2,jj,4,irot),jj=1,7),ch10,&
&       (zarot(3,jj,4,irot),jj=1,7),ch10,&
&       (zarot(4,jj,4,irot),jj=1,7),ch10,&
&       (zarot(5,jj,4,irot),jj=1,7),ch10,&
&       (zarot(6,jj,4,irot),jj=1,7),ch10,&
&       (zarot(7,jj,4,irot),jj=1,7)
       call wrtout(std_out,msg,'COLL')
     end if
   end if

 end do  ! isym loop

end subroutine setsym_ylm
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/setnabla_ylm
!! NAME
!! setnabla_ylm
!!
!! FUNCTION
!! Evaluate several inegrals involving spherical harmonics and their gradient.
!! These integrals are angular part for <phi_i|nabla|phi_j> and <tphi_i|nabla|tphi_j>.
!!
!! INPUTS
!!  mpsang=1+ max. angular momentum
!!
!! OUTPUT
!!  ang_phipphj :: angular part for <phi_i|nabla|phi_j> and <tphi_i|nabla|tphi_j>
!!  ang_phipphj(i,j,1)=\int sin\theta cos\phi Si Sj d\omega
!!  ang_phipphj(i,j,2)=\int cos\theta cos\phi Si \frac{d}{d\theta}Sj d\Omega
!!  ang_phipphj(i,j,3)=\int -sin\phi  Si \frac{d}{d\phi}Sj d\Omega
!!  ang_phipphj(i,j,4)=\int sin\theta sin\phi Si Sj d\Omega
!!  ang_phipphj(i,j,5)=\int cos\theta sin\phi Si \frac{d}{d\theta}Sj d\Omega
!!  ang_phipphj(i,j,6)=\int cos\phi Si \frac{d}{d\phi}Sj d\Omega
!!  ang_phipphj(i,j,7)=\int cos\theta  Si Sj d\Omega
!!  ang_phipphj(i,j,8)=\int -sin\theta Si \frac{d}{d\theta}Sj d\Omega
!!
!!  NOTES
!!   See : Mazevet, S., Torrent, M., Recoules, V. and Jollet, F., High Energy Density Physics, 6, 84-88 (2010)
!!         Calculations of the Transport Properties within the PAW Formalism
!! PARENTS
!!      m_paw_onsite
!!
!! CHILDREN
!!
!! SOURCE

 subroutine setnabla_ylm(ang_phipphj,mpsang)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpsang
!arrays
 real(dp),intent(out) :: ang_phipphj(mpsang**2,mpsang**2,8)

!Local variables-------------------------------
 character(len=500) :: msg
 real(dp) :: ang_phipphj_tmp(16,16,8)

! ************************************************************************

 if (mpsang>4) then
   msg='  Not designed for angular momentum greater than 3 !'
   MSG_ERROR(msg)
 end if

!8 angular integrals for l=0..3, m=-l..+l
!ang_phipphj(1,4,1)=\frac{1}{\sqrt{3}}
!ang_phipphj(2,5,1)=\frac{1}{\sqrt{5}}
!ang_phipphj(3,8,1)=\frac{1}{\sqrt{5}}
!ang_phipphj(4,1,1)=\frac{1}{\sqrt{3}}
!ang_phipphj(4,7,1)=-\frac{1}{\sqrt{15}}
!ang_phipphj(4,9,1)=\frac{1}{\sqrt{5}}
!ang_phipphj(5,2,1)=\frac{1}{\sqrt{5}}
!ang_phipphj(5,10,1)=\sqrt{\frac{3}{14}}
!ang_phipphj(5,12,1)=-\frac{1}{\sqrt{70}}
!ang_phipphj(6,11,1)=\frac{1}{\sqrt{7}}
!ang_phipphj(7,4,1)=-\frac{1}{\sqrt{15}}
!ang_phipphj(7,14,1)=\sqrt{\frac{6}{35}}
!ang_phipphj(8,3,1)=\frac{1}{\sqrt{5}}
!ang_phipphj(8,13,1)=-\sqrt{\frac{3}{35}}
!ang_phipphj(8,15,1)=\frac{1}{\sqrt{7}}
!ang_phipphj(9,4,1)=\frac{1}{\sqrt{5}}
!ang_phipphj(9,14,1)=-\frac{1}{\sqrt{70}}
!ang_phipphj(9,16,1)=\sqrt{\frac{3}{14}}
!ang_phipphj(10,5,1)=\sqrt{\frac{3}{14}}
!ang_phipphj(11,6,1)=\frac{1}{\sqrt{7}}
!ang_phipphj(12,5,1)=-\frac{1}{\sqrt{70}}
!ang_phipphj(13,8,1)=-\sqrt{\frac{3}{35}}
!ang_phipphj(14,7,1)=\sqrt{\frac{6}{35}}
!ang_phipphj(14,9,1)=-\frac{1}{\sqrt{70}}
!ang_phipphj(15,8,1)=\frac{1}{\sqrt{7}}
!ang_phipphj(16,9,1)=\sqrt{\frac{3}{14}}
!ang_phipphj(1,4,2)=\frac{1}{2 \sqrt{3}}
!ang_phipphj(1,14,2)=-\frac{\sqrt{\frac{7}{6}}}{2}
!ang_phipphj(2,5,2)=\frac{1}{2 \sqrt{5}}
!ang_phipphj(3,8,2)=\frac{1}{2 \sqrt{5}}
!ang_phipphj(4,7,2)=-\sqrt{\frac{3}{5}}
!ang_phipphj(4,9,2)=\frac{1}{2 \sqrt{5}}
!ang_phipphj(5,2,2)=\frac{1}{4 \sqrt{5}}
!ang_phipphj(5,10,2)=\frac{\sqrt{\frac{3}{14}}}{2}
!ang_phipphj(5,12,2)=-2 \sqrt{\frac{2}{35}}
!ang_phipphj(6,11,2)=\frac{1}{2 \sqrt{7}}
!ang_phipphj(7,4,2)=\frac{1}{\sqrt{15}}
!ang_phipphj(7,14,2)=\frac{13}{2 \sqrt{210}}
!ang_phipphj(8,3,2)=-\frac{1}{\sqrt{5}}
!ang_phipphj(8,13,2)=-4 \sqrt{\frac{3}{35}}
!ang_phipphj(8,15,2)=\frac{1}{2 \sqrt{7}}
!ang_phipphj(9,4,2)=\frac{1}{4 \sqrt{5}}
!ang_phipphj(9,14,2)=-2 \sqrt{\frac{2}{35}}
!ang_phipphj(9,16,2)=\frac{\sqrt{\frac{3}{14}}}{2}
!ang_phipphj(10,5,2)=\frac{1}{\sqrt{42}}
!ang_phipphj(11,6,2)=-\frac{1}{4 \sqrt{7}}
!ang_phipphj(12,5,2)=\sqrt{\frac{2}{35}}
!ang_phipphj(13,8,2)=2 \sqrt{\frac{3}{35}}
!ang_phipphj(14,7,2)=-2 \sqrt{\frac{6}{35}}
!ang_phipphj(14,9,2)=\sqrt{\frac{2}{35}}
!ang_phipphj(15,8,2)=-\frac{1}{4 \sqrt{7}}
!ang_phipphj(16,9,2)=\frac{1}{\sqrt{42}}
!ang_phipphj(1,4,3)=\frac{\sqrt{3}}{2}
!ang_phipphj(1,14,3)=\frac{\sqrt{\frac{7}{6}}}{2}
!ang_phipphj(2,5,3)=\frac{\sqrt{5}}{2}
!ang_phipphj(3,8,3)=\frac{\sqrt{5}}{2}
!ang_phipphj(4,9,3)=\frac{\sqrt{5}}{2}
!ang_phipphj(5,2,3)=-\frac{\sqrt{5}}{4}
!ang_phipphj(5,10,3)=\frac{\sqrt{\frac{21}{2}}}{2}
!ang_phipphj(6,11,3)=\frac{\sqrt{7}}{2}
!ang_phipphj(7,14,3)=\frac{\sqrt{\frac{35}{6}}}{2}
!ang_phipphj(8,15,3)=\frac{\sqrt{7}}{2}
!ang_phipphj(9,4,3)=-\frac{\sqrt{5}}{4}
!ang_phipphj(9,16,3)=\frac{\sqrt{\frac{21}{2}}}{2}
!ang_phipphj(10,5,3)=-\sqrt{\frac{7}{6}}
!ang_phipphj(11,6,3)=-\frac{\sqrt{7}}{4}
!ang_phipphj(15,8,3)=-\frac{\sqrt{7}}{4}
!ang_phipphj(16,9,3)=-\sqrt{\frac{7}{6}}
!ang_phipphj(1,2,4)=\frac{1}{\sqrt{3}}
!ang_phipphj(2,1,4)=\frac{1}{\sqrt{3}}
!ang_phipphj(2,7,4)=-\frac{1}{\sqrt{15}}
!ang_phipphj(2,9,4)=-\frac{1}{\sqrt{5}}
!ang_phipphj(3,6,4)=\frac{1}{\sqrt{5}}
!ang_phipphj(4,5,4)=\frac{1}{\sqrt{5}}
!ang_phipphj(5,4,4)=\frac{1}{\sqrt{5}}
!ang_phipphj(5,14,4)=-\frac{1}{\sqrt{70}}
!ang_phipphj(5,16,4)=-\sqrt{\frac{3}{14}}
!ang_phipphj(6,3,4)=\frac{1}{\sqrt{5}}
!ang_phipphj(6,13,4)=-\sqrt{\frac{3}{35}}
!ang_phipphj(6,15,4)=-\frac{1}{\sqrt{7}}
!ang_phipphj(7,2,4)=-\frac{1}{\sqrt{15}}
!ang_phipphj(7,12,4)=\sqrt{\frac{6}{35}}
!ang_phipphj(8,11,4)=\frac{1}{\sqrt{7}}
!ang_phipphj(9,2,4)=-\frac{1}{\sqrt{5}}
!ang_phipphj(9,10,4)=\sqrt{\frac{3}{14}}
!ang_phipphj(9,12,4)=\frac{1}{\sqrt{70}}
!ang_phipphj(10,9,4)=\sqrt{\frac{3}{14}}
!ang_phipphj(11,8,4)=\frac{1}{\sqrt{7}}
!ang_phipphj(12,7,4)=\sqrt{\frac{6}{35}}
!ang_phipphj(12,9,4)=\frac{1}{\sqrt{70}}
!ang_phipphj(13,6,4)=-\sqrt{\frac{3}{35}}
!ang_phipphj(14,5,4)=-\frac{1}{\sqrt{70}}
!ang_phipphj(15,6,4)=-\frac{1}{\sqrt{7}}
!ang_phipphj(16,5,4)=-\sqrt{\frac{3}{14}}
!ang_phipphj(1,2,5)=\frac{1}{2 \sqrt{3}}
!ang_phipphj(1,12,5)=-\frac{\sqrt{\frac{7}{6}}}{2}
!ang_phipphj(2,7,5)=-\sqrt{\frac{3}{5}}
!ang_phipphj(2,9,5)=-\frac{1}{2 \sqrt{5}}
!ang_phipphj(3,6,5)=\frac{1}{2 \sqrt{5}}
!ang_phipphj(4,5,5)=\frac{1}{2 \sqrt{5}}
!ang_phipphj(5,4,5)=\frac{1}{4 \sqrt{5}}
!ang_phipphj(5,14,5)=-2 \sqrt{\frac{2}{35}}
!ang_phipphj(5,16,5)=-\frac{\sqrt{\frac{3}{14}}}{2}
!ang_phipphj(6,3,5)=-\frac{1}{\sqrt{5}}
!ang_phipphj(6,13,5)=-4 \sqrt{\frac{3}{35}}
!ang_phipphj(6,15,5)=-\frac{1}{2 \sqrt{7}}
!ang_phipphj(7,2,5)=\frac{1}{\sqrt{15}}
!ang_phipphj(7,12,5)=\frac{13}{2 \sqrt{210}}
!ang_phipphj(8,11,5)=\frac{1}{2 \sqrt{7}}
!ang_phipphj(9,2,5)=-\frac{1}{4 \sqrt{5}}
!ang_phipphj(9,10,5)=\frac{\sqrt{\frac{3}{14}}}{2}
!ang_phipphj(9,12,5)=2 \sqrt{\frac{2}{35}}
!ang_phipphj(10,9,5)=\frac{1}{\sqrt{42}}
!ang_phipphj(11,8,5)=-\frac{1}{4 \sqrt{7}}
!ang_phipphj(12,7,5)=-2 \sqrt{\frac{6}{35}}
!ang_phipphj(12,9,5)=-\sqrt{\frac{2}{35}}
!ang_phipphj(13,6,5)=2 \sqrt{\frac{3}{35}}
!ang_phipphj(14,5,5)=\sqrt{\frac{2}{35}}
!ang_phipphj(15,6,5)=\frac{1}{4 \sqrt{7}}
!ang_phipphj(16,5,5)=-\frac{1}{\sqrt{42}}
!ang_phipphj(1,2,6)=\frac{\sqrt{3}}{2}
!ang_phipphj(1,12,6)=\frac{\sqrt{\frac{7}{6}}}{2}
!ang_phipphj(2,9,6)=-\frac{\sqrt{5}}{2}
!ang_phipphj(3,6,6)=\frac{\sqrt{5}}{2}
!ang_phipphj(4,5,6)=\frac{\sqrt{5}}{2}
!ang_phipphj(5,4,6)=-\frac{\sqrt{5}}{4}
!ang_phipphj(5,16,6)=-\frac{\sqrt{\frac{21}{2}}}{2}
!ang_phipphj(6,15,6)=-\frac{\sqrt{7}}{2}
!ang_phipphj(7,12,6)=\frac{\sqrt{\frac{35}{6}}}{2}
!ang_phipphj(8,11,6)=\frac{\sqrt{7}}{2}
!ang_phipphj(9,2,6)=\frac{\sqrt{5}}{4}
!ang_phipphj(9,10,6)=\frac{\sqrt{\frac{21}{2}}}{2}
!ang_phipphj(10,9,6)=-\sqrt{\frac{7}{6}}
!ang_phipphj(11,8,6)=-\frac{\sqrt{7}}{4}
!ang_phipphj(15,6,6)=\frac{\sqrt{7}}{4}
!ang_phipphj(16,5,6)=\sqrt{\frac{7}{6}}
!ang_phipphj(1,3,7)=\frac{1}{\sqrt{3}}
!ang_phipphj(2,6,7)=\frac{1}{\sqrt{5}}
!ang_phipphj(3,1,7)=\frac{1}{\sqrt{3}}
!ang_phipphj(3,7,7)=\frac{2}{\sqrt{15}}
!ang_phipphj(4,8,7)=\frac{1}{\sqrt{5}}
!ang_phipphj(5,11,7)=\frac{1}{\sqrt{7}}
!ang_phipphj(6,2,7)=\frac{1}{\sqrt{5}}
!ang_phipphj(6,12,7)=2 \sqrt{\frac{2}{35}}
!ang_phipphj(7,3,7)=\frac{2}{\sqrt{15}}
!ang_phipphj(7,13,7)=\frac{3}{\sqrt{35}}
!ang_phipphj(8,4,7)=\frac{1}{\sqrt{5}}
!ang_phipphj(8,14,7)=2 \sqrt{\frac{2}{35}}
!ang_phipphj(9,15,7)=\frac{1}{\sqrt{7}}
!ang_phipphj(11,5,7)=\frac{1}{\sqrt{7}}
!ang_phipphj(12,6,7)=2 \sqrt{\frac{2}{35}}
!ang_phipphj(13,7,7)=\frac{3}{\sqrt{35}}
!ang_phipphj(14,8,7)=2 \sqrt{\frac{2}{35}}
!ang_phipphj(15,9,7)=\frac{1}{\sqrt{7}}
!ang_phipphj(1,3,8)=\frac{2}{\sqrt{3}}
!ang_phipphj(2,6,8)=\frac{3}{\sqrt{5}}
!ang_phipphj(3,7,8)=2 \sqrt{\frac{3}{5}}
!ang_phipphj(4,8,8)=\frac{3}{\sqrt{5}}
!ang_phipphj(5,11,8)=\frac{4}{\sqrt{7}}
!ang_phipphj(6,2,8)=-\frac{1}{\sqrt{5}}
!ang_phipphj(6,12,8)=8 \sqrt{\frac{2}{35}}
!ang_phipphj(7,3,8)=-\frac{2}{\sqrt{15}}
!ang_phipphj(7,13,8)=\frac{12}{\sqrt{35}}
!ang_phipphj(8,4,8)=-\frac{1}{\sqrt{5}}
!ang_phipphj(8,14,8)=8 \sqrt{\frac{2}{35}}
!ang_phipphj(9,15,8)=\frac{4}{\sqrt{7}}
!ang_phipphj(11,5,8)=-\frac{2}{\sqrt{7}}
!ang_phipphj(12,6,8)=-4 \sqrt{\frac{2}{35}}
!ang_phipphj(13,7,8)=-\frac{6}{\sqrt{35}}
!ang_phipphj(14,8,8)=-4 \sqrt{\frac{2}{35}}
!ang_phipphj(15,9,8)=-\frac{2}{\sqrt{7}}


 ang_phipphj_tmp=zero
!
 ang_phipphj_tmp(1,4,1)=0.57735026918962576451_dp
 ang_phipphj_tmp(2,5,1)=0.44721359549995793928_dp
 ang_phipphj_tmp(3,8,1)=0.44721359549995793928_dp
 ang_phipphj_tmp(4,1,1)=0.57735026918962576451_dp
 ang_phipphj_tmp(4,7,1)=-0.25819888974716112568_dp
 ang_phipphj_tmp(4,9,1)=0.44721359549995793928_dp
 ang_phipphj_tmp(5,2,1)=0.44721359549995793928_dp
 ang_phipphj_tmp(5,10,1)=0.46291004988627573078_dp
 ang_phipphj_tmp(5,12,1)=-0.11952286093343936400_dp
 ang_phipphj_tmp(6,11,1)=0.37796447300922722721_dp
 ang_phipphj_tmp(7,4,1)=-0.25819888974716112568_dp
 ang_phipphj_tmp(7,14,1)=0.41403933560541253068_dp
 ang_phipphj_tmp(8,3,1)=0.44721359549995793928_dp
 ang_phipphj_tmp(8,13,1)=-0.29277002188455995381_dp
 ang_phipphj_tmp(8,15,1)=0.37796447300922722721_dp
 ang_phipphj_tmp(9,4,1)=0.44721359549995793928_dp
 ang_phipphj_tmp(9,14,1)=-0.11952286093343936400_dp
 ang_phipphj_tmp(9,16,1)=0.46291004988627573078_dp
 ang_phipphj_tmp(10,5,1)=0.46291004988627573078_dp
 ang_phipphj_tmp(11,6,1)=0.37796447300922722721_dp
 ang_phipphj_tmp(12,5,1)=-0.11952286093343936400_dp
 ang_phipphj_tmp(13,8,1)=-0.29277002188455995381_dp
 ang_phipphj_tmp(14,7,1)=0.41403933560541253068_dp
 ang_phipphj_tmp(14,9,1)=-0.11952286093343936400_dp
 ang_phipphj_tmp(15,8,1)=0.37796447300922722721_dp
 ang_phipphj_tmp(16,9,1)=0.46291004988627573078_dp
!
 ang_phipphj_tmp(1,4,2)=0.28867513459481288225_dp
 ang_phipphj_tmp(1,14,2)=-0.54006172486732168591_dp
 ang_phipphj_tmp(2,5,2)=0.22360679774997896964_dp
 ang_phipphj_tmp(3,8,2)=0.22360679774997896964_dp
 ang_phipphj_tmp(4,7,2)=-0.77459666924148337704_dp
 ang_phipphj_tmp(4,9,2)=0.22360679774997896964_dp
 ang_phipphj_tmp(5,2,2)=0.11180339887498948482_dp
 ang_phipphj_tmp(5,10,2)=0.23145502494313786539_dp
 ang_phipphj_tmp(5,12,2)=-0.47809144373375745599_dp
 ang_phipphj_tmp(6,11,2)=0.18898223650461361361_dp
 ang_phipphj_tmp(7,4,2)=0.25819888974716112568_dp
 ang_phipphj_tmp(7,14,2)=0.44854261357253024157_dp
 ang_phipphj_tmp(8,3,2)=-0.44721359549995793928_dp
 ang_phipphj_tmp(8,13,2)=-1.1710800875382398152_dp
 ang_phipphj_tmp(8,15,2)=0.18898223650461361361_dp
 ang_phipphj_tmp(9,4,2)=0.11180339887498948482_dp
 ang_phipphj_tmp(9,14,2)=-0.47809144373375745599_dp
 ang_phipphj_tmp(9,16,2)=0.23145502494313786539_dp
 ang_phipphj_tmp(10,5,2)=0.15430334996209191026_dp
 ang_phipphj_tmp(11,6,2)=-0.094491118252306806804_dp
 ang_phipphj_tmp(12,5,2)=0.23904572186687872799_dp
 ang_phipphj_tmp(13,8,2)=0.58554004376911990761_dp
 ang_phipphj_tmp(14,7,2)=-0.82807867121082506136_dp
 ang_phipphj_tmp(14,9,2)=0.23904572186687872799_dp
 ang_phipphj_tmp(15,8,2)=-0.094491118252306806804_dp
 ang_phipphj_tmp(16,9,2)=0.15430334996209191026_dp
!
 ang_phipphj_tmp(1,4,3)=0.86602540378443864676_dp
 ang_phipphj_tmp(1,14,3)=0.54006172486732168591_dp
 ang_phipphj_tmp(2,5,3)=1.1180339887498948482_dp
 ang_phipphj_tmp(3,8,3)=1.1180339887498948482_dp
 ang_phipphj_tmp(4,9,3)=1.1180339887498948482_dp
 ang_phipphj_tmp(5,2,3)=-0.55901699437494742410_dp
 ang_phipphj_tmp(5,10,3)=1.6201851746019650577_dp
 ang_phipphj_tmp(6,11,3)=1.3228756555322952953_dp
 ang_phipphj_tmp(7,14,3)=1.2076147288491198811_dp
 ang_phipphj_tmp(8,15,3)=1.3228756555322952953_dp
 ang_phipphj_tmp(9,4,3)=-0.55901699437494742410_dp
 ang_phipphj_tmp(9,16,3)=1.6201851746019650577_dp
 ang_phipphj_tmp(10,5,3)=-1.0801234497346433718_dp
 ang_phipphj_tmp(11,6,3)=-0.66143782776614764763_dp
 ang_phipphj_tmp(15,8,3)=-0.66143782776614764763_dp
 ang_phipphj_tmp(16,9,3)=-1.0801234497346433718_dp
!
 ang_phipphj_tmp(1,2,4)=0.57735026918962576451_dp
 ang_phipphj_tmp(2,1,4)=0.57735026918962576451_dp
 ang_phipphj_tmp(2,7,4)=-0.25819888974716112568_dp
 ang_phipphj_tmp(2,9,4)=-0.44721359549995793928_dp
 ang_phipphj_tmp(3,6,4)=0.44721359549995793928_dp
 ang_phipphj_tmp(4,5,4)=0.44721359549995793928_dp
 ang_phipphj_tmp(5,4,4)=0.44721359549995793928_dp
 ang_phipphj_tmp(5,14,4)=-0.11952286093343936400_dp
 ang_phipphj_tmp(5,16,4)=-0.46291004988627573078_dp
 ang_phipphj_tmp(6,3,4)=0.44721359549995793928_dp
 ang_phipphj_tmp(6,13,4)=-0.29277002188455995381_dp
 ang_phipphj_tmp(6,15,4)=-0.37796447300922722721_dp
 ang_phipphj_tmp(7,2,4)=-0.25819888974716112568_dp
 ang_phipphj_tmp(7,12,4)=0.41403933560541253068_dp
 ang_phipphj_tmp(8,11,4)=0.37796447300922722721_dp
 ang_phipphj_tmp(9,2,4)=-0.44721359549995793928_dp
 ang_phipphj_tmp(9,10,4)=0.46291004988627573078_dp
 ang_phipphj_tmp(9,12,4)=0.11952286093343936400_dp
 ang_phipphj_tmp(10,9,4)=0.46291004988627573078_dp
 ang_phipphj_tmp(11,8,4)=0.37796447300922722721_dp
 ang_phipphj_tmp(12,7,4)=0.41403933560541253068_dp
 ang_phipphj_tmp(12,9,4)=0.11952286093343936400_dp
 ang_phipphj_tmp(13,6,4)=-0.29277002188455995381_dp
 ang_phipphj_tmp(14,5,4)=-0.11952286093343936400_dp
 ang_phipphj_tmp(15,6,4)=-0.37796447300922722721_dp
 ang_phipphj_tmp(16,5,4)=-0.46291004988627573078_dp
!
 ang_phipphj_tmp(1,2,5)=0.28867513459481288225_dp
 ang_phipphj_tmp(1,12,5)=-0.54006172486732168591_dp
 ang_phipphj_tmp(2,7,5)=-0.77459666924148337704_dp
 ang_phipphj_tmp(2,9,5)=-0.22360679774997896964_dp
 ang_phipphj_tmp(3,6,5)=0.22360679774997896964_dp
 ang_phipphj_tmp(4,5,5)=0.22360679774997896964_dp
 ang_phipphj_tmp(5,4,5)=0.11180339887498948482_dp
 ang_phipphj_tmp(5,14,5)=-0.47809144373375745599_dp
 ang_phipphj_tmp(5,16,5)=-0.23145502494313786539_dp
 ang_phipphj_tmp(6,3,5)=-0.44721359549995793928_dp
 ang_phipphj_tmp(6,13,5)=-1.1710800875382398152_dp
 ang_phipphj_tmp(6,15,5)=-0.18898223650461361361_dp
 ang_phipphj_tmp(7,2,5)=0.25819888974716112568_dp
 ang_phipphj_tmp(7,12,5)=0.44854261357253024157_dp
 ang_phipphj_tmp(8,11,5)=0.18898223650461361361_dp
 ang_phipphj_tmp(9,2,5)=-0.11180339887498948482_dp
 ang_phipphj_tmp(9,10,5)=0.23145502494313786539_dp
 ang_phipphj_tmp(9,12,5)=0.47809144373375745599_dp
 ang_phipphj_tmp(10,9,5)=0.15430334996209191026_dp
 ang_phipphj_tmp(11,8,5)=-0.094491118252306806804_dp
 ang_phipphj_tmp(12,7,5)=-0.82807867121082506136_dp
 ang_phipphj_tmp(12,9,5)=-0.23904572186687872799_dp
 ang_phipphj_tmp(13,6,5)=0.58554004376911990761_dp
 ang_phipphj_tmp(14,5,5)=0.23904572186687872799_dp
 ang_phipphj_tmp(15,6,5)=0.094491118252306806804_dp
 ang_phipphj_tmp(16,5,5)=-0.15430334996209191026_dp
!
 ang_phipphj_tmp(1,2,6)=0.86602540378443864676_dp
 ang_phipphj_tmp(1,12,6)=0.54006172486732168591_dp
 ang_phipphj_tmp(2,9,6)=-1.1180339887498948482_dp
 ang_phipphj_tmp(3,6,6)=1.1180339887498948482_dp
 ang_phipphj_tmp(4,5,6)=1.1180339887498948482_dp
 ang_phipphj_tmp(5,4,6)=-0.55901699437494742410_dp
 ang_phipphj_tmp(5,16,6)=-1.6201851746019650577_dp
 ang_phipphj_tmp(6,15,6)=-1.3228756555322952953_dp
 ang_phipphj_tmp(7,12,6)=1.2076147288491198811_dp
 ang_phipphj_tmp(8,11,6)=1.3228756555322952953_dp
 ang_phipphj_tmp(9,2,6)=0.55901699437494742410_dp
 ang_phipphj_tmp(9,10,6)=1.6201851746019650577_dp
 ang_phipphj_tmp(10,9,6)=-1.0801234497346433718_dp
 ang_phipphj_tmp(11,8,6)=-0.66143782776614764763_dp
 ang_phipphj_tmp(15,6,6)=0.66143782776614764763_dp
 ang_phipphj_tmp(16,5,6)=1.0801234497346433718_dp
!
 ang_phipphj_tmp(1,3,7)=0.57735026918962576451_dp
 ang_phipphj_tmp(2,6,7)=0.44721359549995793928_dp
 ang_phipphj_tmp(3,1,7)=0.57735026918962576451_dp
 ang_phipphj_tmp(3,7,7)=0.51639777949432225136_dp
 ang_phipphj_tmp(4,8,7)=0.44721359549995793928_dp
 ang_phipphj_tmp(5,11,7)=0.37796447300922722721_dp
 ang_phipphj_tmp(6,2,7)=0.44721359549995793928_dp
 ang_phipphj_tmp(6,12,7)=0.47809144373375745599_dp
 ang_phipphj_tmp(7,3,7)=0.51639777949432225136_dp
 ang_phipphj_tmp(7,13,7)=0.50709255283710994651_dp
 ang_phipphj_tmp(8,4,7)=0.44721359549995793928_dp
 ang_phipphj_tmp(8,14,7)=0.47809144373375745599_dp
 ang_phipphj_tmp(9,15,7)=0.37796447300922722721_dp
 ang_phipphj_tmp(11,5,7)=0.37796447300922722721_dp
 ang_phipphj_tmp(12,6,7)=0.47809144373375745599_dp
 ang_phipphj_tmp(13,7,7)=0.50709255283710994651_dp
 ang_phipphj_tmp(14,8,7)=0.47809144373375745599_dp
 ang_phipphj_tmp(15,9,7)=0.37796447300922722721_dp
!
 ang_phipphj_tmp(1,3,8)=1.1547005383792515290_dp
 ang_phipphj_tmp(2,6,8)=1.3416407864998738178_dp
 ang_phipphj_tmp(3,7,8)=1.5491933384829667541_dp
 ang_phipphj_tmp(4,8,8)=1.3416407864998738178_dp
 ang_phipphj_tmp(5,11,8)=1.5118578920369089089_dp
 ang_phipphj_tmp(6,2,8)=-0.44721359549995793928_dp
 ang_phipphj_tmp(6,12,8)=1.9123657749350298240_dp
 ang_phipphj_tmp(7,3,8)=-0.51639777949432225136_dp
 ang_phipphj_tmp(7,13,8)=2.0283702113484397860_dp
 ang_phipphj_tmp(8,4,8)=-0.44721359549995793928_dp
 ang_phipphj_tmp(8,14,8)=1.9123657749350298240_dp
 ang_phipphj_tmp(9,15,8)=1.5118578920369089089_dp
 ang_phipphj_tmp(11,5,8)=-0.75592894601845445443_dp
 ang_phipphj_tmp(12,6,8)=-0.95618288746751491198_dp
 ang_phipphj_tmp(13,7,8)=-1.0141851056742198930_dp
 ang_phipphj_tmp(14,8,8)=-0.95618288746751491198_dp
 ang_phipphj_tmp(15,9,8)=-0.75592894601845445443_dp

 ang_phipphj(:,:,:)=ang_phipphj_tmp(1:mpsang**2,1:mpsang**2,:)

 end subroutine setnabla_ylm
!!***
!----------------------------------------------------------------------

!!****f* m_paw_sphharm/rfactorial
!! NAME
!! rfactorial
!!
!! FUNCTION
!! Private function
!! Calculates N! as a double precision real.
!!
!! INPUTS
!!   nn=number to use
!!
!! OUTPUT
!!   factorial= n! (real)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental function rfactorial(nn)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: nn
 real(dp) :: rfactorial

!Local variables ---------------------------------------
!scalars
 integer :: ii

! *********************************************************************

 rfactorial=one
 do ii=2,nn
   rfactorial=rfactorial*ii
 end do

end function rfactorial
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/perms
!! NAME
!! perms
!!
!! FUNCTION
!! Private function
!! Returns N!/(N-k)!  if N>=0 and N>k ; otherwise 0 is returned
!!
!! INPUTS
!!   kk=number k to use
!!   nn=number N to use
!!
!! OUTPUT
!!   perms= n!/(n-k)!
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function perms(nn,kk)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: kk,nn
 real(dp) :: perms

!Local variables ---------------------------------------
!scalars
 integer :: ii
 real(dp) :: pp

! *********************************************************************

 if (nn>=0.and.nn>=kk) then
   pp=1._dp
   do ii=nn-kk+1,nn
     pp=pp*ii
   end do
 else
   pp=0._dp
 end if

 perms=pp

end function perms
!!***

!----------------------------------------------------------------------
!!****f* m_pawang/gaunt
!! NAME
!! gaunt
!!
!! FUNCTION
!! Returns gaunt coefficient, i.e.
!!   the integral of Sqrt[4 \pi] Y*(l_i,m_i) Y*(ll,mm) Y(l_j,m_j)
!!   See the 3-j and 6-j symbols by Rotenberg, etc., (Technology Press, 1959), pg.5.
!!
!! INPUTS
!!   ll,mm,l1,l2,m1,m2= six quantum numbers defining the Gaunt coef.
!!
!! OUTPUT
!!   gaunt(ll,mm,l1,l2,m1,m2)=the value of the integral
!!
!! CHILDREN
!!
!! SOURCE

function gaunt(ll,mm,l1,m1,l2,m2)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: l1,l2,ll,m1,m2,mm
 real(dp) :: gaunt

!Local variables ------------------------------
!scalars
 integer :: i1,i2,j1,j1half,j2,j2half,j3,j3half,j_half,jj,k1,k2,n1,n2
 real(dp) :: argument,sign,sum,xx,yy
 logical :: ok

!************************************************************************

 gaunt=zero;sum=zero;ok =.true.

 if((-m1-mm+m2) /= 0) ok = .false.
 if(abs(m1) > l1) ok = .false.
 if(abs(mm) > ll) ok = .false.
 if(abs(m2) > l2) ok = .false.

 jj = l1 + ll + l2
 if (mod(jj,2)/=0) ok = .false.
 j1 = jj-2*l2
 j2 = jj-2*ll
 j3 = jj-2*l1

 if (j1<0 .or. j2<0 .or. j3<0) ok = .false.

 if (ok) then

   xx = (2 * l1 + 1) * (2 * ll + 1) * (2 * l2 + 1)

   j1half = j1/2
   j2half = j2/2
   j3half = j3/2
   j_half = jj/2

   gaunt = (-1)**j1half * sqrt(xx)
   gaunt = gaunt * rfactorial(j2)*rfactorial(j3)/rfactorial(jj+1)
   gaunt = gaunt * rfactorial(j_half)/(rfactorial(j1half)&
&                * rfactorial(j2half)*rfactorial(j3half))

   yy = rfactorial(l2 + m2) * rfactorial(l2 - m2)

   if (mm>=0) then
     yy = yy * perms(ll+mm,2*mm)
   else
     yy = yy / perms(ll-mm,-2*mm)
   end if

   if (m1>=0) then
     yy = yy / perms(l1+m1,2*m1)
   else
     yy = yy * perms(l1-m1,-2*m1)
   end if

   gaunt = gaunt * sqrt(yy)

   i1 = l2 - ll - m1
   i2 = l2 - l1 + mm
   k1 = -min(0, i1, i2)
   n1 = l1 + m1
   n2 = ll - mm
   k2 = min(j1, n1, n2)

   sign = 1._dp
   if(k1>0) sign = (-1._dp)**k1

   argument = sign     * perms(n1,k1)/rfactorial(k1)
   argument = argument * perms(n2,k1)/rfactorial(i1 + k1)
   argument = argument * perms(j1,k1)/rfactorial(i2 + k1)
   sum = sum + argument

   sign = -sign
   k1 = k1 + 1
   do while(k1 <= k2)
     argument = sign     * perms(n1, k1)/rfactorial(k1)
     argument = argument * perms(n2, k1)/rfactorial(i1 + k1)
     argument = argument * perms(j1, k1)/rfactorial(i2 + k1)
     sum = sum + argument
     sign = -sign
     k1 = k1 + 1
   end do

 end if

 gaunt = gaunt * sum

 end function gaunt
!!***

!----------------------------------------------------------------------

!!****f* m_pawang/realgaunt
!! NAME
!! realgaunt
!!
!! FUNCTION
!! This routine compute "real Gaunt coefficients", i.e. gaunt
!! coefficients according to "real spherical harmonics"
!!
!! INPUTS
!!  l_max= max. value of ang. momentum l+1;  Gaunt coeffs up to
!!          [(2*l_max-1,m),(l_max,m),(l_max,m)] are computed
!!
!! OUTPUT
!!  gntselect((2*l_max-1)**2,l_max**2*(l_max**2+1)/2)=
!!          selection rules for Gaunt coefficients
!!          if Gaunt coeff. is zero, gntselect=0
!!          if Gaunt coeff. is non-zero, gntselect is the index of
!!                           the coeff. in realgnt(:) array
!!  ngnt= number of non-zero Gaunt coefficients
!!  realgnt((2*l_max-1)**2*l_max**4)= non-zero real Gaunt coefficients
!!
!! PARENTS
!!      m_paw_slater,m_pawang,m_pawpwij
!!
!! CHILDREN
!!
!! SOURCE

subroutine realgaunt(l_max,ngnt,gntselect,realgnt)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: l_max
 integer,intent(out) :: ngnt
!arrays
 integer,intent(out) :: gntselect((2*l_max-1)**2,l_max**2*(l_max**2+1)/2)
 real(dp),intent(out) :: realgnt((2*l_max-1)**2*(l_max)**4)

!Local variables ------------------------------
!scalars
 integer :: ilm1,ilm2,ilmp1,k0lm1,klm1,l1,l2,ll,lp1,m1,m2,mm,mm1,mm2,mm3,mp1
 real(dp) :: c11,c12,c21,c22,c31,c32,fact,realgnt_tmp
!arrays
 integer,allocatable :: ssgn(:)
 type(coeff3_type), allocatable :: coeff(:)

!************************************************************************

!Initialize output arrays with zeros.
 gntselect = 0; realgnt = zero

!Compute matrix cc where Sl=cc*Yl (Sl=real sph. harm.)
!------------------------------------------------
 LIBPAW_DATATYPE_ALLOCATE(coeff,(4*l_max-3))
 do ll=1,4*l_max-3
   LIBPAW_ALLOCATE(coeff(ll)%value,(2,2*ll-1,2*ll-1))
   coeff(ll)%value(:,:,:)=zero
   coeff(ll)%value(1,ll,ll)=one
   do mm=1,ll-1
     coeff(ll)%value(1,ll+mm,ll+mm)= (-1._dp)**mm/sqrt(2._dp)
     coeff(ll)%value(1,ll-mm,ll+mm)= ( 1._dp)    /sqrt(2._dp)
     coeff(ll)%value(2,ll+mm,ll-mm)=-(-1._dp)**mm/sqrt(2._dp)
     coeff(ll)%value(2,ll-mm,ll-mm)= ( 1._dp)    /sqrt(2._dp)
   end do
 end do

 LIBPAW_ALLOCATE(ssgn,(l_max**2))
 ssgn(:)=1
 if (l_max>0) then
   do l1=1,l_max-1
     ilm1=1+l1**2+l1
     do m1=-l1,-1
       ssgn(ilm1+m1)=-1
     end do
   end do
 end if

 ngnt=0

!Loop on (lp1,mp1)
!------------------------------------------------
 do lp1=0,l_max-1
   do mp1=-lp1,lp1
     ilmp1=1+lp1**2+lp1+mp1
     k0lm1=ilmp1*(ilmp1-1)/2

!    Loop on (l1,m1)<=(lp1,mp1)
!    ------------------------------------------------
     do l1=0,l_max-1
       do m1=-l1,l1
         ilm1=1+l1**2+l1+m1

         if (ilm1<=ilmp1) then

           klm1=k0lm1+ilm1
           gntselect(:,klm1)=0

!          Loop on (l2,m2)
!          ------------------------------------------------
           do l2=abs(l1-lp1),l1+lp1,2
             do m2=-l2,l2
               ilm2=1+l2**2+l2+m2

!              Real Gaunt coeffs selection rules
!              ------------------------------------------------
               if ((l2<=l1+lp1).and.&
&               (((m1== mp1).and.((m2==0).or.(m2==2*abs(mp1)))).or.&
&               ((m1==-mp1).and.(m2==-abs(m1)-abs(mp1))).or.&
&               ((abs(m1)/=(abs(mp1)).and.&
&               ((m2==ssgn(ilm1)*ssgn(ilmp1)*   (abs(m1)+abs(mp1))).or.&
&               (m2==ssgn(ilm1)*ssgn(ilmp1)*abs(abs(m1)-abs(mp1)))&
               ))))) then

!                Compute selected real Gaunt coefficient
!                ------------------------------------------------
                 realgnt_tmp=zero
                 do mm1=-l1,l1
                   c11=coeff(l1+1)%value(1,l1+mm1+1,l1+m1+1)
                   c12=coeff(l1+1)%value(2,l1+mm1+1,l1+m1+1)
                   do mm2= -lp1,lp1
                     c21=coeff(lp1+1)%value(1,lp1+mm2+1,lp1+mp1+1)
                     c22=coeff(lp1+1)%value(2,lp1+mm2+1,lp1+mp1+1)
                     do mm3= -l2,l2
                       c31=coeff(l2+1)%value(1,l2+mm3+1,l2+m2+1)
                       c32=coeff(l2+1)%value(2,l2+mm3+1,l2+m2+1)
                       fact=c11*c21*c31  -  c12*c22*c31&
&                       -c11*c22*c32  -  c12*c21*c32
                       if((abs(fact)>=tol12).and.(mm3==-mm2-mm1)) &
&                       realgnt_tmp=realgnt_tmp+fact*(-1)**mm2 &
&                       *gaunt(l2,mm3,l1,mm1,lp1,-mm2)
                     end do
                   end do
                 end do

!                Count additional non-zero real Gaunt coeffs
!                ------------------------------------------------
                 if (abs(realgnt_tmp)>=tol12) then
                   ngnt=ngnt+1
                   gntselect(ilm2,klm1)=ngnt
                   realgnt(ngnt)=realgnt_tmp/sqrt(four_pi)
                 end if

!                End loops
!                ------------------------------------------------
               end if
             end do
           end do
         end if
       end do
     end do
   end do
 end do

!Deallocate memory
!------------------------------------------------
 do ll=1,4*l_max-3
   LIBPAW_DEALLOCATE(coeff(ll)%value)
 end do
 LIBPAW_DATATYPE_DEALLOCATE(coeff)
 LIBPAW_DEALLOCATE(ssgn)

end subroutine realgaunt
!!***

!----------------------------------------------------------------------

!!****f* m_paw_sphharm/nablarealgaunt
!! NAME
!! nablarealgaunt
!!
!! FUNCTION
!! Evaluate inegrals involving spherical harmonics and their gradient.
!! These integrals are angular part for <nablaphi|nablaphj> and <tnablaphi|tnablaphj>.
!!
!! INPUTS
!!  mpsang=1+ max. angular momentum
!!
!! OUTPUT
!! nnablagnt= number of non zero integrals
!! nabgauntselect(mpsang**2,mpsang**2,mpsang**2)= stores the index of the non zero integrals
!! nablagaunt((mpsang**2)**3)= stores the integrals' values
!!  NOTES
!!  
!! PARENTS
!!      m_pawang
!!
!! CHILDREN
!!
!! SOURCE

subroutine nablarealgaunt(mpsang,nnablagnt,nabgauntselect,nablagaunt) 
 implicit none

 !Arguments ---------------------------------------------
 !scalars
 integer, intent(in) :: mpsang
 integer, intent(out) :: nnablagnt
 !array
 integer,intent(out) :: nabgauntselect(mpsang**2,mpsang**2,mpsang**2)
 real(dp),intent(out) :: nablagaunt((mpsang**2)**3)

 !Local variables ---------------------------------------
 !array


nabgauntselect(:,:,:)=-1
nablagaunt(:)=zero



if (mpsang>=2) then
  nabgauntselect(1,2,2)=1  ; nablagaunt(1)=0.5641895835477563_dp !(1/sqrt(pi)) 
  nabgauntselect(1,3,3)=2  ; nablagaunt(2)=0.5641895835477563_dp !(1/sqrt(pi))
  nabgauntselect(1,4,4)=3  ; nablagaunt(3)=0.5641895835477563_dp !(1/sqrt(pi))
  nnablagnt=3
end if

if (mpsang>=3) then


 nabgauntselect(7,2,2)=4  ; nablagaunt(4)=0.126156626101008_dp !\frac{1}{2\sqrt{5\pi}}
 nabgauntselect(9,2,2)=5 ; nablagaunt(5)=0.2185096861184158_dp !\frac{0.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(6,2,3)=6 ; nablagaunt(6)=-0.2185096861184158_dp !-\frac{0.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(5,2,4)=7 ; nablagaunt(7)=-0.2185096861184158_dp !-\frac{0.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(4,2,5)=8 ; nablagaunt(8)=0.6555290583552474_dp !\frac{1.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(3,2,6)=9 ; nablagaunt(9)=0.6555290583552474_dp !\frac{1.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(2,2,7)=10 ; nablagaunt(10)=-0.37846987830302403_dp !\frac{-3}{2\sqrt{5\pi}}
 nabgauntselect(2,2,9)=11 ; nablagaunt(11)=-0.6555290583552474_dp !\frac{-1.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(6,3,2)=12 ; nablagaunt(12)=-0.2185096861184158_dp !-\frac{0.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(7,3,3)=13 ; nablagaunt(13)=-0.252313252202016_dp !-\frac{1}{\sqrt{5\pi}}
 nabgauntselect(8,3,4)=14 ; nablagaunt(14)=-0.2185096861184158_dp !-\frac{0.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(2,3,6)=15 ; nablagaunt(15)=0.6555290583552474_dp !\frac{1.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(3,3,7)=16 ; nablagaunt(16)=0.75693974607408354_dp !\frac{3}{\sqrt{5\pi}
 nabgauntselect(4,3,8)=17 ; nablagaunt(17)=0.6555290583552474_dp !\frac{1.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(5,4,2)=18 ; nablagaunt(18)=-0.2185096861184158_dp !-\frac{0.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(8,4,3)=19 ; nablagaunt(19)=-0.2185096861184158_dp !-\frac{0.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(7,4,4)=20 ; nablagaunt(20)=0.126156626101008_dp !\frac{1}{2\sqrt{5\pi}}
 nabgauntselect(9,4,4)=21 ; nablagaunt(21)=-0.2185096861184158_dp !-\frac{0.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(2,4,5)=22 ; nablagaunt(22)=0.6555290583552474_dp !\frac{1.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(4,4,7)=23 ; nablagaunt(23)=-0.37846987830302403_dp !\frac{-3}{2\sqrt{5\pi}}
 nabgauntselect(3,4,8)=24 ; nablagaunt(24)=0.6555290583552474_dp !\frac{1.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(4,4,9)=25 ; nablagaunt(25)=0.6555290583552474_dp !\frac{1.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(4,5,2)=26 ; nablagaunt(26)=0.6555290583552474_dp !\frac{1.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(2,5,4)=27 ; nablagaunt(27)=0.6555290583552474_dp !\frac{1.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(1,5,5)=28 ; nablagaunt(28)=1.692568750643269_dp !\frac{3}{\sqrt{\pi}}
 nabgauntselect(7,5,5)=29 ; nablagaunt(29)=-0.5406712547186058_dp !-\frac{3}{7}\sqrt{\frac{5}{\pi}}
 nabgauntselect(8,5,6)=30 ; nablagaunt(30)=0.4682350416823196_dp !\frac{3}{14}\sqrt{\frac{15}{\pi}}
 nabgauntselect(5,5,7)=31 ; nablagaunt(31)=-0.5406712547186058_dp !-\frac{3}{7}\sqrt{\frac{5}{\pi}}
 nabgauntselect(6,5,8)=32 ; nablagaunt(32)=0.4682350416823196_dp !\frac{3}{14}\sqrt{\frac{15}{\pi}}
 nabgauntselect(3,6,2)=33 ; nablagaunt(33)=0.6555290583552474_dp !\frac{1.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(2,6,3)=34 ; nablagaunt(34)=0.6555290583552474_dp !\frac{1.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(8,6,5)=35 ; nablagaunt(35)=0.4682350416823196_dp !\frac{3}{14}\sqrt{\frac{15}{\pi}}
 nabgauntselect(1,6,6)=36 ; nablagaunt(36)=1.692568750643269_dp !\dfrac{3}{\sqrt{\pi}}
 nabgauntselect(7,6,6)=37 ; nablagaunt(37)=0.2703356273593029_dp !\frac{3}{14}\sqrt{\frac{5}{\pi}}
 nabgauntselect(9,6,6)=38 ; nablagaunt(38)=-0.4682350416823196_dp !\frac{3}{14}\sqrt{\frac{15}{\pi}}
 nabgauntselect(6,6,7)=39 ; nablagaunt(39)=0.2703356273593029_dp !\frac{3}{14}\sqrt{\frac{5}{\pi}}
 nabgauntselect(5,6,8)=40 ; nablagaunt(40)=0.4682350416823196_dp !\frac{3}{14}\sqrt{\frac{15}{\pi}}
 nabgauntselect(6,6,9)=41 ; nablagaunt(41)=-0.4682350416823196_dp !-\frac{3}{14}\sqrt{\frac{15}{\pi}}
 nabgauntselect(2,7,2)=42 ; nablagaunt(42)=-0.37846987830302403_dp!\frac{-3}{2\sqrt{5\pi}}
 nabgauntselect(3,7,3)=43 ; nablagaunt(43)=0.7569397566060481_dp !\frac{3}{\sqrt{5\pi}}
 nabgauntselect(4,7,4)=44 ; nablagaunt(44)=-0.37846987830302403_dp!\frac{-3}{2\sqrt{5\pi}}
 nabgauntselect(5,7,5)=45 ; nablagaunt(45)=-0.5406712547186058_dp !-\frac{3}{7}\sqrt{\frac{5}{\pi}}
 nabgauntselect(6,7,6)=46 ; nablagaunt(46)=0.2703356273593029_dp !\frac{3}{14}\sqrt{\frac{5}{\pi}}
 nabgauntselect(1,7,7)=47 ; nablagaunt(47)=1.692568750643269_dp !\frac{3}{\sqrt{\pi}}
 nabgauntselect(7,7,7)=48 ; nablagaunt(48)=0.5406712547186058_dp !\frac{3}{7}\sqrt{\frac{5}{\pi}}
 nabgauntselect(8,7,8)=49 ; nablagaunt(49)=0.2703356273593029_dp !\frac{3}{14}\sqrt{\frac{5}{\pi}}
 nabgauntselect(9,7,9)=50 ; nablagaunt(50)=-0.5406712547186058_dp !-\frac{3}{7}\sqrt{\frac{5}{\pi}}
 nabgauntselect(4,8,3)=51 ; nablagaunt(51)=0.6555290583552474_dp  !\frac{3}{2}\sqrt{\frac{3}{5\pi}}
 nabgauntselect(3,8,4)=52 ; nablagaunt(52)=0.6555290583552474_dp !\frac{3}{2}\sqrt{\frac{3}{5\pi}}
 nabgauntselect(6,8,5)=53 ; nablagaunt(53)=0.4682350416823196_dp !\frac{3}{14}\sqrt{\frac{15}{\pi}}
 nabgauntselect(5,8,6)=54 ; nablagaunt(54)=0.4682350416823196_dp !\frac{3}{14}\sqrt{\frac{15}{\pi}}
 nabgauntselect(8,8,7)=55 ; nablagaunt(55)=0.2703356273593029_dp !\frac{3}{14}\sqrt{\frac{5}{\pi}}
 nabgauntselect(1,8,8)=56 ; nablagaunt(56)=1.692568750643269_dp !\frac{3}{\sqrt{\pi}}
 nabgauntselect(7,8,8)=57 ; nablagaunt(57)=0.2703356273593029_dp !\frac{3}{14}\sqrt{\frac{5}{\pi}}
 nabgauntselect(9,8,8)=58 ; nablagaunt(58)=0.4682350416823196_dp !\frac{3}{14}\sqrt{\frac{15}{\pi}}
 nabgauntselect(8,8,9)=59 ; nablagaunt(59)=0.4682350416823196_dp !\frac{3}{14}\sqrt{\frac{15}{\pi}}
 nabgauntselect(9,2,2)=60 ; nablagaunt(60)=0.2185096861184158_dp !!!!! \frac{0.5\sqrt{3}}{\sqrt{5\pi}}
 nabgauntselect(9,4,4)=61 ; nablagaunt(61)=0.6555290583552474_dp !\frac{3}{2}\sqrt{\frac{3}{5\pi}}
 nabgauntselect(6,9,6)=62 ; nablagaunt(62)=-0.4682350416823196_dp !-\frac{3}{14}\sqrt{\frac{15}{\pi}}
 nabgauntselect(9,9,7)=63 ; nablagaunt(63)=-0.5406712547186058_dp !-\frac{3}{7}\sqrt{\frac{5}{\pi}}
 nabgauntselect(8,9,8)=64 ; nablagaunt(64)=0.4682350416823196_dp !\frac{3}{14}\sqrt{\frac{15}{\pi}}
 nabgauntselect(1,9,9)=65 ; nablagaunt(65)=1.692568750643269_dp !\frac{3}{\sqrt{\pi}}
 nabgauntselect(7,9,9)=66 ; nablagaunt(66)=-0.5406712547186058_dp !-\frac{3}{7}\sqrt{\frac{5}{\pi}}
 nabgauntselect(4,9,4)=67 ; nablagaunt(67)=0.6555290583552474_dp !\frac{3}{2}\sqrt{\frac{3}{5\pi}}
 nnablagnt=67
 end if

end subroutine nablarealgaunt


END MODULE m_paw_sphharm
!!***

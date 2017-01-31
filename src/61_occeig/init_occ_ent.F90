!{\src2tex{textfont=tt}}
!!****f* ABINIT/init_occ_ent
!! NAME
!! init_occ_ent
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2016 ABINIT group (the_author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! PARENTS
!!      getnel
!!
!! CHILDREN
!!      spline
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine init_occ_ent(entfun,limit,nptsdiv2,occfun,occopt,option,smdfun,tphysel,tsmear,tsmearinv,xgrid)

 use defs_basis
 use m_splines
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_occ_ent'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: occopt,option
 real(dp),intent(in) :: tphysel,tsmear
 integer,intent(inout) :: nptsdiv2
 real(dp),intent(out) :: limit,tsmearinv
 real(dp),intent(inout) :: entfun(-nptsdiv2:nptsdiv2,2),occfun(-nptsdiv2:nptsdiv2,2)
 real(dp),intent(inout) :: smdfun(-nptsdiv2:nptsdiv2,2),xgrid(-nptsdiv2:nptsdiv2)


!Local variables-------------------------------
!scalars
 integer :: algo,ii,jj,nconvd2
 integer :: nmaxFD,nminFD
 integer,parameter :: nptsdiv2_def=6000
 integer,save :: dblsmr,occopt_prev=-9999
 real(dp),parameter :: maxFDarg=500.0_dp
 real(dp),save :: convlim,incconv,limit_occ,tphysel_prev=-9999,tsmear_prev=-9999
 real(dp) :: aa,dsqrpi,encorr,factor
 real(dp) :: expinc,expx22,expxo2,gauss,increm
 real(dp) :: resFD1,resFD2,resFD3,resFD4,resmom,resmom1,resmom2
 real(dp) :: resmom3,resmom4,secmom,smom1,smom2,thdmom,tmom1,tmom2,tmpexpsum
 real(dp) :: tmpsmdfun,tratio,tt,xx,yp1,ypn
 character(len=500) :: message
!arrays
 real(dp),save :: entfun_prev(-nptsdiv2_def:nptsdiv2_def,2),occfun_prev(-nptsdiv2_def:nptsdiv2_def,2)
 real(dp),save :: smdfun_prev(-nptsdiv2_def:nptsdiv2_def,2),xgrid_prev(-nptsdiv2_def:nptsdiv2_def)
 real(dp),allocatable :: entder(:),occder(:),smd1(:),smd2(:)
 real(dp),allocatable :: smdder(:),tgrid(:),work(:),workfun(:)

! *************************************************************************

!Initialize the occupation function and generalized entropy function,
!at the beginning, or if occopt changed

 if(option==-1)then
   nptsdiv2 = nptsdiv2_def
   return
 end if


 if(occopt_prev/=occopt           .or. &
& abs(tsmear_prev-tsmear)  >tol12 .or. &
& abs(tphysel_prev-tphysel)>tol12       ) then
!  write(std_out,*) 'INIT_OCC_ENT CHANGE ..........'
   occopt_prev=occopt
   tsmear_prev=tsmear
   tphysel_prev=tphysel

!  Check whether input values of tphysel tsmear and occopt are consistent
   dblsmr = 0
   if (abs(tphysel)>tol12) then
!    Use re-smearing scheme
     if (abs(tsmear)>tol12) then
       dblsmr = 1
!      Use FD occupations (one smearing) only with "physical" temperature tphysel
     else if (occopt /= 3) then
       write(message, '(a,i6,a)' )' tphysel /= 0, tsmear == 0, but occopt is not = 3, but ',occopt,'.'
       MSG_ERROR(message)
     end if
   end if
!  write(std_out,*) 'getnel : input read.'
!  write(std_out,*) '  dblsmr = ', dblsmr
!  write(std_out,*) '  tphysel, tsmear = ', tphysel, tsmear

   ABI_ALLOCATE(entder,(-nptsdiv2_def:nptsdiv2_def))
   ABI_ALLOCATE(occder,(-nptsdiv2_def:nptsdiv2_def))
   ABI_ALLOCATE(smdder,(-nptsdiv2_def:nptsdiv2_def))
   ABI_ALLOCATE(workfun,(-nptsdiv2_def:nptsdiv2_def))
   ABI_ALLOCATE(work,(-nptsdiv2_def:nptsdiv2_def))

!  Prepare the points on the grid
!  limit is the value of the argument that will give 0.0 or 1.0 , with
!  less than about 1.0d-15 error for 4<=occopt<=8, and less than about 1.0d-12
!  error for occopt==3. It is not worth to compute the function beyond
!  that point. Even with a less severe requirement, it is significantly
!  larger for occopt==3, with an exponential
!  tail, than for the other occupation functions, with a Gaussian tail.
!  Note that these values are useful in newocc.f also.
   limit_occ=6.0_dp
   if(occopt==3)limit_occ=30.0_dp
   if(dblsmr /= 0) then
     tratio = tsmear / tphysel
     limit_occ=30.0_dp + 6.0_dp*tratio
   end if

!  With nptsdiv2_def=6000 (thus increm=0.001 for 4<=occopt<=8,
!  and increm=0.005 for occopt==3, the O(1/N4) algorithm gives 1.0d-12
!  accuracy on the stored values occfun and entfun. These, together
!  with smdfun and xgrid_prev, need permanently about 0.67 MB, which is affordable.
   increm=limit_occ/nptsdiv2_def
   do ii=-nptsdiv2_def,nptsdiv2_def
     xgrid_prev(ii)=ii*increm
   end do

!  ---------------------------------------------------------
!  Ordinary (unique) smearing function
!  ---------------------------------------------------------
   if (dblsmr == 0) then

!    Compute the unnormalized smeared delta function between -limit_occ and +limit_occ
!    (well, they are actually normalized ...)
     if(occopt==3)then

!      Fermi-Dirac
       do ii=0,nptsdiv2_def
         xx=xgrid_prev(ii)
         smdfun_prev( ii,1)=0.25_dp/(cosh(xx/2.0_dp)**2)
         smdfun_prev(-ii,1)=smdfun_prev(ii,1)
       end do

     else if(occopt==4 .or. occopt==5)then

!      Cold smearing of Marzari, two values of the "a" parameter being possible
!      first value gives minimization of the bump
       if(occopt==4)aa=-.5634
!      second value gives monotonic occupation function
       if(occopt==5)aa=-.8165

       dsqrpi=1.0_dp/sqrt(pi)
       do ii=0,nptsdiv2_def
         xx=xgrid_prev(ii)
         gauss=dsqrpi*exp(-xx**2)
         smdfun_prev( ii,1)=gauss*(1.5_dp+xx*(-aa*1.5_dp+xx*(-1.0_dp+aa*xx)))
         smdfun_prev(-ii,1)=gauss*(1.5_dp+xx*( aa*1.5_dp+xx*(-1.0_dp-aa*xx)))
       end do

     else if(occopt==6)then

!      First order Hermite-Gaussian of Paxton and Methfessel
       dsqrpi=1.0_dp/sqrt(pi)
       do ii=0,nptsdiv2_def
         xx=xgrid_prev(ii)
         smdfun_prev( ii,1)=dsqrpi*(1.5_dp-xx**2)*exp(-xx**2)
         smdfun_prev(-ii,1)=smdfun_prev(ii,1)
       end do

     else if(occopt==7)then

!      Gaussian smearing
       dsqrpi=1.0_dp/sqrt(pi)
       do ii=0,nptsdiv2_def
         xx=xgrid_prev(ii)
         smdfun_prev( ii,1)=dsqrpi*exp(-xx**2)
         smdfun_prev(-ii,1)=smdfun_prev(ii,1)
       end do

     else if(occopt==8)then

!      Constant value of the delta function over the smearing interval, for testing purposes only.
       do ii=0,nptsdiv2_def
         xx=xgrid_prev(ii)
         if(xx>half+tol8)then
           smdfun_prev( ii,1)=zero
         else if(xx<half-tol8)then
           smdfun_prev( ii,1)=one
         else
           smdfun_prev( ii,1)=half
         end if
         smdfun_prev(-ii,1)=smdfun_prev(ii,1)
       end do

     else
       write(message, '(a,i0,a)' )' Occopt=',occopt,' is not allowed in getnel. '
       MSG_BUG(message)
     end if

!    ---------------------------------------------------------
!    smear FD delta with occopt delta calculated in smdfun_prev
!    ---------------------------------------------------------
   else if (dblsmr /= 0) then

     nconvd2 = 6000
     convlim = 10.0_dp
     incconv = convlim / nconvd2

!    store smearing functions in smd1 and smd2
     ABI_ALLOCATE(smd1,(-nconvd2:nconvd2))
     ABI_ALLOCATE(smd2,(-nconvd2:nconvd2))
     ABI_ALLOCATE(tgrid,(-nconvd2:nconvd2))

!    FD function in smd1( ii) and second smearing delta in smd2( ii)
!
!    smd1(:) contains delta_FD ( x )
     do ii=0,nconvd2
       tgrid(ii)=ii*incconv
       tgrid(-ii)=-tgrid(ii)
       tt=tgrid(ii)
       smd1( ii)=0.25_dp/(cosh(tt/2.0_dp)**2)
       smd1(-ii)=smd1(ii)
     end do

!    check input values of occopt and fill smd2(:) with appropriate data:
!    smd2(:) contains delta_resmear ( x )
     if(occopt == 3) then
       write(message, '(a,a)' )&
&       'Occopt=3 is not allowed as a re-smearing.', &
&       'Use a single FD, or re-smear with a different delta type (faster cutoff). '
       MSG_ERROR(message)
     else if(occopt==4 .or. occopt==5)then
!      Cold smearing of Marzari, two values of the "a" parameter being possible
!      first value gives minimization of the bump
       if(occopt==4)aa=-.5634
!      second value gives monotonic occupation function
       if(occopt==5)aa=-.8165

       dsqrpi=1.0_dp/sqrt(pi)
       do ii=0,nconvd2
         tt=tgrid(ii)
         gauss=dsqrpi*exp(-tt**2)
         smd2( ii)=gauss*(1.5_dp+tt*(-aa*1.5_dp+tt*(-1.0_dp+aa*tt)))
         smd2(-ii)=gauss*(1.5_dp+tt*( aa*1.5_dp+tt*(-1.0_dp-aa*tt)))
       end do
     else if(occopt==6)then
       dsqrpi=1.0_dp/sqrt(pi)
       do ii=0,nconvd2
         tt=tgrid(ii)
         smd2( ii)=dsqrpi*(1.5_dp-tt**2)*exp(-tt**2)
         smd2(-ii)=smd2(ii)
       end do
     else if(occopt==7)then
       dsqrpi=1.0_dp/sqrt(pi)
       do ii=0,nconvd2
         tt=tgrid(ii)
         smd2( ii)=dsqrpi*exp(-tt**2)
         smd2(-ii)=smd2(ii)
       end do
     else if(occopt==8)then
       do ii=0,nconvd2
         tt=tgrid(ii)
         if(tt>half+tol8)then
           smd2( ii)=zero
         else if(tt<half-tol8)then
           smd2( ii)=one
         else
           smd2( ii)=half
         end if
         smd2(-ii)=smd2(ii)
       end do
     else
       write(message, '(a,i0,a)' )' Occopt= ',occopt,' is not allowed in getnel. '
       MSG_BUG(message)
     end if


!    Use O(1/N4) algorithm from Num Rec (see below)
!
!    The grid for the convoluted delta is taken (conservatively)
!    to be that for the FD delta ie 6000 pts in [-limit_occ;limit_occ]
!    Smearing functions are given on [-dbllim;dbllim] and the grid must
!    superpose the normal grid on [-limit_occ:limit_occ]
!    The maximal interval for integration of the convolution is
!    [-dbllim+limit_occ+lim(delta2);dbllim-limit_occ-lim(delta2)] =
!    [-dbllim+36;dbllim-36]

!    test the smdFD function for extreme values:
!    do jj=-nptsdiv2_def,-nptsdiv2_def
!    do ii=-nconvd2+4,nconvd2
!    call smdFD(xgrid_prev(jj) - tgrid(ii)*tratio, resFD)
!    write(std_out,*) 'ii jj = ', ii,jj, ' smdFD (', &
!    &    xgrid_prev(jj) - tgrid(ii)*tratio, ') ', resFD
!    end do
!    end do

     expinc = exp(half*incconv*tratio)

!    jj = position of point at which we are calculating smdfun_prev
     do jj=-nptsdiv2_def,nptsdiv2_def
!      Do not care about the 8 boundary points,
!      where the values should be extremely small anyway
       smdfun_prev(jj,1)=0.0_dp
!      only add contribution with delta_FD > 1.0d-100
       nmaxFD = floor  (( maxFDarg+xgrid_prev(jj)) / tratio / incconv )
       nmaxFD = min (nmaxFD, nconvd2)
       nminFD = ceiling((-maxFDarg+xgrid_prev(jj)) / tratio / incconv )
       nminFD = max (nminFD, -nconvd2)

!      Calculate the Fermi-Dirac distrib at point xgrid_prev(jj)-tgrid(ii)*tratio
       expxo2 = exp (-half*(xgrid_prev(jj) - (nminFD)*incconv*tratio))
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD4 = tmpexpsum * tmpexpsum
       expxo2 = expxo2*expinc
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD3 = tmpexpsum * tmpexpsum
       expxo2 = expxo2*expinc
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD2 = tmpexpsum * tmpexpsum
       expxo2 = expxo2*expinc
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD1 = tmpexpsum * tmpexpsum

!      core contribution to the integral with constant weight (48)
       tmpsmdfun = 0.0_dp
       do ii=nminFD+4,nmaxFD-4
         expxo2 = expxo2*expinc
!        tmpexpsum = 1.0_dp / (expxo2 + 1.0_dp / expxo2 )
         expx22 = expxo2*expxo2
         tmpexpsum = expxo2 / (expx22 + 1.0_dp)
         tmpsmdfun = tmpsmdfun + smd2(ii) * tmpexpsum * tmpexpsum
       end do

!      Add on end contributions for show (both functions smd and smdFD are very small
       smdfun_prev(jj,1)=smdfun_prev(jj,1)       +48.0_dp*tmpsmdfun             &
&       + 31.0_dp*smd2(nminFD+3)*resFD1 -11.0_dp*smd2(nminFD+2)*resFD2 &
&       +  5.0_dp*smd2(nminFD+1)*resFD3 -       smd2(nminFD)*resFD4

       expxo2 = expxo2*expinc
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD1 = tmpexpsum * tmpexpsum
       expxo2 = expxo2*expinc
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD2 = tmpexpsum * tmpexpsum
       expxo2 = expxo2*expinc
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD3 = tmpexpsum * tmpexpsum
       expxo2 = expxo2*expinc
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD4 = tmpexpsum * tmpexpsum

!      Contribution above
       smdfun_prev(jj,1)=smdfun_prev(jj,1)                                      &
&       + 31.0_dp*smd2(nmaxFD-3)*resFD1  -11.0_dp*smd2(nmaxFD-2)*resFD2 &
&       +  5.0_dp*smd2(nmaxFD-1)*resFD3  -       smd2(nmaxFD)*resFD4
       smdfun_prev(jj,1)=incconv*smdfun_prev(jj,1)/48.0_dp
     end do

     secmom = 0.0_dp
     thdmom = 0.0_dp
     resmom4 = xgrid_prev(-nptsdiv2_def  )*xgrid_prev(-nptsdiv2_def  )*smdfun_prev(-nptsdiv2_def  ,  1)
     resmom3 = xgrid_prev(-nptsdiv2_def+1)*xgrid_prev(-nptsdiv2_def+1)*smdfun_prev(-nptsdiv2_def+1,  1)
     resmom2 = xgrid_prev(-nptsdiv2_def+2)*xgrid_prev(-nptsdiv2_def+2)*smdfun_prev(-nptsdiv2_def+2,  1)
     resmom1 = xgrid_prev(-nptsdiv2_def+3)*xgrid_prev(-nptsdiv2_def+3)*smdfun_prev(-nptsdiv2_def+3,  1)
     resmom  = xgrid_prev(-nptsdiv2_def+4)*xgrid_prev(-nptsdiv2_def+4)*smdfun_prev(-nptsdiv2_def+4,  1)
     do ii=-nptsdiv2_def+4,nptsdiv2_def-1
       secmom = secmom +                                   &
&       ( 17.0_dp*xgrid_prev(ii)  *xgrid_prev(ii)  *smdfun_prev(ii,  1)   &
&       +42.0_dp*xgrid_prev(ii-1)*xgrid_prev(ii-1)*smdfun_prev(ii-1,1)   &
&       -16.0_dp*xgrid_prev(ii-2)*xgrid_prev(ii-2)*smdfun_prev(ii-2,1)   &
&       + 6.0_dp*xgrid_prev(ii-3)*xgrid_prev(ii-3)*smdfun_prev(ii-3,1)   &
&       -       xgrid_prev(ii-4)*xgrid_prev(ii-4)*smdfun_prev(ii-4,1)  )
       resmom4 = resmom3
       resmom3 = resmom2
       resmom2 = resmom1
       resmom1 = resmom
       resmom  = xgrid_prev(ii+1)  *xgrid_prev(ii+1)  *smdfun_prev(ii+1,  1)
     end do
     secmom=increm * secmom / 48.0_dp
!    thdmom=increm * thdmom / 48.0_dp
!
!    smom1  = second moment of delta in smd1(:)
!    smom2  = second moment of delta in smd2(:)
!
     smom1  = 0.0_dp
     smom2  = 0.0_dp
     tmom1  = 0.0_dp
     tmom2  = 0.0_dp
     do ii=-nconvd2+4,nconvd2
       smom1 = smom1+                                       &
&       ( 17.0_dp*tgrid(ii)  *tgrid(ii)  *smd1(ii)         &
&       +42.0_dp*tgrid(ii-1)*tgrid(ii-1)*smd1(ii-1)       &
&       -16.0_dp*tgrid(ii-2)*tgrid(ii-2)*smd1(ii-2)       &
&       + 6.0_dp*tgrid(ii-3)*tgrid(ii-3)*smd1(ii-3)       &
&       -       tgrid(ii-4)*tgrid(ii-4)*smd1(ii-4)  )
       smom2 = smom2+                                       &
&       ( 17.0_dp*tgrid(ii)  *tgrid(ii)  *smd2(ii  )     &
&       +42.0_dp*tgrid(ii-1)*tgrid(ii-1)*smd2(ii-1)     &
&       -16.0_dp*tgrid(ii-2)*tgrid(ii-2)*smd2(ii-2)     &
&       + 6.0_dp*tgrid(ii-3)*tgrid(ii-3)*smd2(ii-3)     &
&       -       tgrid(ii-4)*tgrid(ii-4)*smd2(ii-4)  )
     end do
     smom1 =incconv * smom1  / 48.0_dp
     smom2 =incconv * smom2  / 48.0_dp
!    tmom1 =incconv * tmom1  / 48.0_dp
!    tmom2 =incconv * tmom2  / 48.0_dp

     encorr =  smom2*tratio*tratio/secmom

!    DEBUG
!    write(std_out,*) ' getnel : debug, secmoms = ', secmom, smom1, smom2
!    write(std_out,*) ' getnel : debug, thdmoms = ', thdmom, tmom1, tmom2
!    write(std_out,*) ' getnel : encorr = ', encorr
!    ENDDEBUG

     ABI_DEALLOCATE(tgrid)
     ABI_DEALLOCATE(smd1)
     ABI_DEALLOCATE(smd2)

   end if

!  --------------------------------------------------------
!  end of smearing function initialisation, dblsmr case
!  --------------------------------------------------------


!  Now that the smeared delta function has been initialized, compute the
!  occupation function
   occfun_prev(-nptsdiv2_def,1)=zero
   entfun_prev(-nptsdiv2_def,1)=zero

!  Different algorithms are possible, corresponding to the formulas
!  (4.1.11), (4.1.12) and (4.1.14) in Numerical recipes (pp 107 and 108),
!  with respective O(1/N2), O(1/N3), O(1/N4) convergence, where N is the
!  number of points in the interval.
   algo=4

   if(algo==2)then

!    Extended trapezoidal rule (4.1.11), taken in a cumulative way
     do ii=-nptsdiv2_def+1,nptsdiv2_def
       occfun_prev(ii,1)=occfun_prev(ii-1,1)+increm*(smdfun_prev(ii,1)+smdfun_prev(ii-1,1))/2.0_dp
       entfun_prev(ii,1)=entfun_prev(ii-1,1)+increm*&
&       ( -xgrid_prev(ii)*smdfun_prev(ii,1) -xgrid_prev(ii-1)*smdfun_prev(ii-1,1) )/2.0_dp
     end do

   else if(algo==3)then

!    Derived from (4.1.12). Converges as O(1/N3).
!    Do not care about the following points,
!    where the values are extremely small anyway
     occfun_prev(-nptsdiv2_def+1,1)=0.0_dp ;   entfun_prev(-nptsdiv2_def+1,1)=0.0_dp
     do ii=-nptsdiv2_def+2,nptsdiv2_def
       occfun_prev(ii,1)=occfun_prev(ii-1,1)+increm*&
&       ( 5.0_dp*smdfun_prev(ii,1) + 8.0_dp*smdfun_prev(ii-1,1) - smdfun_prev(ii-2,1) )/12.0_dp
       entfun_prev(ii,1)=entfun_prev(ii-1,1)+increm*&
&       ( 5.0_dp*(-xgrid_prev(ii)  )*smdfun_prev(ii,1)  &
&       +8.0_dp*(-xgrid_prev(ii-1))*smdfun_prev(ii-1,1)&
&       -      (-xgrid_prev(ii-2))*smdfun_prev(ii-2,1) )/12.0_dp
     end do

   else if(algo==4)then

!    Derived from (4.1.14)- alternative extended Simpsons rule. Converges as O(1/N4).
!    Do not care about the following points,
!    where the values are extremely small anyway
     occfun_prev(-nptsdiv2_def+1,1)=0.0_dp ;   entfun_prev(-nptsdiv2_def+1,1)=0.0_dp
     occfun_prev(-nptsdiv2_def+2,1)=0.0_dp ;   entfun_prev(-nptsdiv2_def+2,1)=0.0_dp
     occfun_prev(-nptsdiv2_def+3,1)=0.0_dp ;   entfun_prev(-nptsdiv2_def+3,1)=0.0_dp
     do ii=-nptsdiv2_def+4,nptsdiv2_def
       occfun_prev(ii,1)=occfun_prev(ii-1,1)+increm*&
&       ( 17.0_dp*smdfun_prev(ii,1)  &
&       +42.0_dp*smdfun_prev(ii-1,1)&
&       -16.0_dp*smdfun_prev(ii-2,1)&
&       + 6.0_dp*smdfun_prev(ii-3,1)&
&       -       smdfun_prev(ii-4,1) )/48.0_dp
       entfun_prev(ii,1)=entfun_prev(ii-1,1)+increm*&
&       ( 17.0_dp*(-xgrid_prev(ii)  )*smdfun_prev(ii,1)  &
&       +42.0_dp*(-xgrid_prev(ii-1))*smdfun_prev(ii-1,1)&
&       -16.0_dp*(-xgrid_prev(ii-2))*smdfun_prev(ii-2,1)&
&       + 6.0_dp*(-xgrid_prev(ii-3))*smdfun_prev(ii-3,1)&
&       -       (-xgrid_prev(ii-4))*smdfun_prev(ii-4,1) )/48.0_dp
     end do

   end if ! End of choice between different algorithms for integration

!  Normalize the functions (actually not needed for occopt=3..7)
   factor=1.0_dp/occfun_prev(nptsdiv2_def,1)
   smdfun_prev(:,1)=smdfun_prev(:,1)*factor
   occfun_prev(:,1)=occfun_prev(:,1)*factor
   entfun_prev(:,1)=entfun_prev(:,1)*factor

!  Compute the cubic spline fitting of the smeared delta function
   yp1=0.0_dp ; ypn=0.0_dp
   workfun(:)=smdfun_prev(:,1)
   call spline(xgrid_prev, workfun, (2*nptsdiv2_def+1), yp1, ypn, smdder)
   smdfun_prev(:,2)=smdder(:)

!  Compute the cubic spline fitting of the occupation function
   yp1=0.0_dp ; ypn=0.0_dp
   workfun(:)=occfun_prev(:,1)
   call spline(xgrid_prev, workfun, (2*nptsdiv2_def+1), yp1, ypn, occder)
   occfun_prev(:,2)=occder(:)

!  Compute the cubic spline fitting of the entropy function
   yp1=0.0_dp ; ypn=0.0_dp
   workfun(:)=entfun_prev(:,1)
   call spline(xgrid_prev, workfun, (2*nptsdiv2_def+1), yp1, ypn, entder)
   entfun_prev(:,2)=entder(:)

   ABI_DEALLOCATE(entder)
   ABI_DEALLOCATE(occder)
   ABI_DEALLOCATE(smdder)
   ABI_DEALLOCATE(work)
   ABI_DEALLOCATE(workfun)

 end if

 if (abs(tphysel)<tol12) then
   tsmearinv=one/tsmear
 else
   tsmearinv=one/tphysel
 end if

 entfun(:,:) = entfun_prev(:,:)
 occfun(:,:) = occfun_prev(:,:)
 smdfun(:,:) = smdfun_prev(:,:)
 xgrid(:) = xgrid_prev(:)
 limit = limit_occ
 nptsdiv2 = nptsdiv2_def

end subroutine init_occ_ent
!!***

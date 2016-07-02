!{\src2tex{textfont=tt}}
!!****f* ABINIT/vgh_rho
!! NAME
!! vgh_rho
!!
!! FUNCTION
!! The general procedure to obtain the value, the gradient and the hessian
!! of the density of electrons in the point vv (in cart.coord).
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! INPUTS
!! vv(3)=position
!! chs  1 only valence density
!!      2 only core density
!!      0 total density
!!     -2 iat, ipos are nulify and ignored
!!     -1 iat,ipos = index of atom if vv is inside
!!         the "core sphere (rminl)", 0 otherwise
!!
!! OUTPUT
!! rho,grho(3),hrho(3,3) - density, gradient of density, hessian of density
!!                                 (cart. coord)
!! iat, ipos - index of the nearest atom (except chs < 0 see above)
!! rdmin  - the distance to the nearest atom
!!
!! SIDE EFFECTS
!!  This routine also works on the data contained in the defs_aimprom and defs_aimfields modules
!!
!! PARENTS
!!      addout,aim_follow,cpdrv,critic,critics,integrho,onestep,plint,rsurf
!!
!! CHILDREN
!!      bschg1,bschg2
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine vgh_rho(vv,rho,grho,hrho,rdmin,iat,ipos,chs)

 use m_profiling_abi

 use defs_basis
 use defs_parameters
 use defs_aimprom
 use defs_aimfields

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vgh_rho'
 use interfaces_63_bader, except_this_one => vgh_rho
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chs
 integer,intent(inout) :: iat,ipos
 real(dp),intent(out) :: rdmin,rho
!arrays
 real(dp),intent(in) :: vv(3)
 real(dp),intent(out) :: grho(3),hrho(3,3)

!Local variables ------------------------------
!scalars
 integer :: ii,inmax,inmin,inx,jj,kk,ll,nn,oii,omm,onn
 integer :: selct
! real(dp),save :: cumul_cpu=0.0_dp,cumul_cpu_old=0.0_dp
 real(dp),save :: tcpui,tcpuo,twalli
 real(dp),save :: twallo
 real(dp) :: aa,bb,cc,cgrad1_rr_inv,coeff,dd,rr,rr2,rr_inv
 real(dp) :: rrad2_nn,rrad_nn,ss,uu,uu_inv,val,vt1,vt2,vt3,vw1,vw2
! real(dp) :: ss_inv
 real(dp) :: vw3
!arrays
 integer :: indx(3),inii(4,3)
 real(dp) :: cgrad(3),ches(3,3),cof(2,3),ddstar(6),ddu(2),grd(4)
 real(dp) :: hh(4,2),hrh(2),lder(4),pom2sq(2,3),pomsq(2,3)
 real(dp) :: rhstar(6),sqder(6,4),sqvlr(6,4),trsf(3,3),xx(3)
 real(dp),pointer :: ptddx(:,:,:),ptddy(:,:,:),ptddz(:,:,:),ptrho(:,:,:)

!************************************************************************
 tcpui=0.0_dp
 tcpuo=0.0_dp
 twalli=0.0_dp
 twallo=0.0_dp

 nullify(ptddx,ptddy,ptddz,ptrho)

 selct=chs

 if (selct/=2) then

!  call timein(tcpui,twalli)

!  TRANSFORMATION TO THE REDUCED COORD.

   xx(:)=vv(:)
   call bschg1(xx,-1)

!  call timein(tcpuo,twallo)
!  cumul_cpu=cumul_cpu+(tcpuo-tcpui)

!  REDUCTION TO THE PRIMITIVE CELL

   do ii=1,3
     if (xx(ii) >= one-tol12 ) then
       xx(ii)=xx(ii)-aint(xx(ii))
     elseif (xx(ii) < -tol12 ) then
       xx(ii)=xx(ii)-floor(xx(ii))
     end if
   end do


!  DETERMINATION OF THE INDEX IN THE GRID

   do ii=1,3
     indx(ii)=aint(xx(ii)*ngfft(ii))
     bb=(xx(ii)-indx(ii)*dix(ii))*ngfft(ii)
     if (indx(ii)==ngfft(ii)) then
       indx(ii)=1
       xx(ii)=0._dp
     else
       indx(ii)=indx(ii)+1
     end if

!    Explicit handling to avoid numeric problems

     if (bb > 1._dp+tol12 ) then
       cof(1,ii)=0._dp
       cof(2,ii)=1._dp
     elseif (bb < -tol12 ) then
       cof(1,ii)=1._dp
       cof(2,ii)=0._dp
     else
       cof(1,ii)=1._dp-bb
       cof(2,ii)=bb
     end if
   end do

!  3D INTERPOLATION OF THE VALENCE DENSITY

!  determination of the values of density and of its second derivative
!  at the "star" = constructed at vv with primitive directions
!  To interpolation the values at the faces of the grid cell are needed

   rhstar(:)=0._dp
   sqder(:,:)=0._dp
   sqvlr(:,:)=0._dp
   ddstar(:)=0._dp
   pomsq(:,:)=0._dp
   pom2sq(:,:)=0._dp

   oii=1; onn=1; omm=1
   if (indx(1)==ngfft(1)) oii=1-ngfft(1)
   if (indx(2)==ngfft(2)) onn=1-ngfft(2)
   if (indx(3)==ngfft(3)) omm=1-ngfft(3)

!  the values in the corners of the grid cell

   ptddx=>ddx(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm)
   ptddy=>ddy(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm)
   ptddz=>ddz(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm)
   ptrho=>dvl(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm)

!  the coefficients for spline interpolation of density and its derivation
   do ii=1,3
     do jj=1,2
       pomsq(jj,ii)=(cof(jj,ii)*cof(jj,ii)*cof(jj,ii)-cof(jj,ii))/6._dp*dix(ii)*dix(ii)
       pom2sq(jj,ii)=(3._dp*cof(jj,ii)*cof(jj,ii)-1._dp)/6._dp*dix(ii)
       if (jj==1) pom2sq(jj,ii)=-pom2sq(jj,ii)
     end do
   end do


   do ii=1,2
     do jj=1,2
       do kk=1,2
         ddstar(ii)=ddstar(ii)+cof(jj,2)*cof(kk,3)*ptddx(ii,jj,kk)
         ddstar(ii+2)=ddstar(ii+2)+cof(jj,3)*cof(kk,1)*ptddy(kk,ii,jj)
         ddstar(ii+4)=ddstar(ii+4)+cof(jj,1)*cof(kk,2)*ptddz(jj,kk,ii)
         sqder(ii,jj)=sqder(ii,jj)+cof(kk,2)*ptddz(ii,kk,jj)
         sqder(ii,jj+2)=sqder(ii,jj+2)+cof(kk,3)*ptddy(ii,jj,kk)
         sqder(ii+2,jj)=sqder(ii+2,jj)+cof(kk,3)*ptddx(jj,ii,kk)
         sqder(ii+2,jj+2)=sqder(ii+2,jj+2)+cof(kk,1)*ptddz(kk,ii,jj)
         sqder(ii+4,jj)=sqder(ii+4,jj)+cof(kk,1)*ptddy(kk,jj,ii)
         sqder(ii+4,jj+2)=sqder(ii+4,jj+2)+cof(kk,2)*ptddx(jj,kk,ii)
         sqvlr(ii,jj)=sqvlr(ii,jj)+cof(kk,2)*ptrho(ii,kk,jj)+pomsq(kk,2)*ptddy(ii,kk,jj)
         sqvlr(ii,jj+2)=sqvlr(ii,jj+2)+cof(kk,3)*ptrho(ii,jj,kk)+pomsq(kk,3)*ptddz(ii,jj,kk)
         sqvlr(ii+2,jj+2)=sqvlr(ii+2,jj+2)+cof(kk,1)*ptrho(kk,ii,jj)+pomsq(kk,1)*ptddx(kk,ii,jj)
       end do
     end do
   end do

   do ii=1,2
     do jj=1,2
       sqvlr(ii+2,jj)=sqvlr(jj,ii+2)
       sqvlr(ii+4,jj)=sqvlr(jj+2,ii+2)
       sqvlr(ii+4,jj+2)=sqvlr(jj,ii)
     end do
   end do

   do ii=1,2
     do jj=1,2
       rhstar(ii)=rhstar(ii)+cof(jj,3)*sqvlr(ii,jj)+pomsq(jj,3)*sqder(ii,jj)+&
&       cof(jj,2)*sqvlr(ii,jj+2)+pomsq(jj,2)*sqder(ii,jj+2)
       rhstar(ii+2)=rhstar(ii+2)+cof(jj,1)*sqvlr(ii+2,jj)+pomsq(jj,1)*sqder(ii+2,jj)+&
&       cof(jj,3)*sqvlr(ii+2,jj+2)+pomsq(jj,3)*sqder(ii+2,jj+2)
       rhstar(ii+4)=rhstar(ii+4)+cof(jj,2)*sqvlr(ii+4,jj)+pomsq(jj,2)*sqder(ii+4,jj)+&
&       cof(jj,1)*sqvlr(ii+4,jj+2)+pomsq(jj,1)*sqder(ii+4,jj+2)
     end do
   end do
   rhstar(:)=rhstar(:)/2._dp

   rho=0._dp
   grho(:)=0._dp
   hrho(:,:)=0._dp
   kk=1; nn=1
   do ii=1,5,2
     do jj=1,2
       nn=-nn
       rho=rho+cof(jj,kk)*rhstar(ii+jj-1)+pomsq(jj,kk)*ddstar(ii+jj-1)
       grho(kk)=grho(kk)+pom2sq(jj,kk)*ddstar(ii+jj-1)
       hrho(kk,kk)=hrho(kk,kk)+cof(jj,kk)*ddstar(ii+jj-1)
       grho(kk)=grho(kk)+nn*rhstar(ii+jj-1)/dix(kk)
     end do
     kk=kk+1
   end do
   rho=rho/3._dp

!  Off-diagonal elements of the hessian

!  for the speed reasons the polynomial interpolation
!  for second derivation fields is used in this case
!  but the last step is always done by spline interpolation.


   do ii=1,3
     do jj=-1,2
       inii(jj+2,ii)=indx(ii)+jj
       if (inii(jj+2,ii) < 1) inii(jj+2,ii)=inii(jj+2,ii)+ngfft(ii)
       if (inii(jj+2,ii) > ngfft(ii)) inii(jj+2,ii)=inii(jj+2,ii)-ngfft(ii)
     end do
   end do

!  Not very nice

   do ii=1,3
     select case (ii)
     case (1)
       do jj=1,4
         ddu(1)=cof(1,2)*ddz(inii(jj,1),inii(2,2),inii(2,3))+cof(2,2)*ddz(inii(jj,1),inii(3,2),inii(2,3))
         ddu(2)=cof(1,2)*ddz(inii(jj,1),inii(2,2),inii(3,3))+cof(2,2)*ddz(inii(jj,1),inii(3,2),inii(3,3))
         hrh(1)=cof(1,2)*dvl(inii(jj,1),inii(2,2),inii(2,3))+cof(2,2)*dvl(inii(jj,1),inii(3,2),inii(2,3))+&
&         pomsq(1,2)*ddy(inii(jj,1),inii(2,2),inii(2,3))+pomsq(2,2)*ddy(inii(jj,1),inii(3,2),inii(2,3))
         hrh(2)=cof(1,2)*dvl(inii(jj,1),inii(2,2),inii(3,3))+cof(2,2)*dvl(inii(jj,1),inii(3,2),inii(3,3))+&
&         pomsq(1,2)*ddy(inii(jj,1),inii(2,2),inii(3,3))+pomsq(2,2)*ddy(inii(jj,1),inii(3,2),inii(3,3))
         hh(jj,2)=(hrh(2)-hrh(1))/dix(3)+pom2sq(1,3)*ddu(1)+pom2sq(2,3)*ddu(2)

         ddu(1)=cof(1,3)*ddy(inii(jj,1),inii(2,2),inii(2,3))+cof(2,3)*ddy(inii(jj,1),inii(2,2),inii(3,3))
         ddu(2)=cof(1,3)*ddy(inii(jj,1),inii(3,2),inii(2,3))+cof(2,3)*ddy(inii(jj,1),inii(3,2),inii(3,3))
         hrh(1)=cof(1,3)*dvl(inii(jj,1),inii(2,2),inii(2,3))+cof(2,3)*dvl(inii(jj,1),inii(2,2),inii(3,3))+&
&         pomsq(1,3)*ddz(inii(jj,1),inii(2,2),inii(2,3))+pomsq(2,3)*ddz(inii(jj,1),inii(2,2),inii(3,3))
         hrh(2)=cof(1,3)*dvl(inii(jj,1),inii(3,2),inii(2,3))+cof(2,3)*dvl(inii(jj,1),inii(3,2),inii(3,3))+&
&         pomsq(1,3)*ddz(inii(jj,1),inii(3,2),inii(2,3))+pomsq(2,3)*ddz(inii(jj,1),inii(3,2),inii(3,3))
         hh(jj,1)=(hrh(2)-hrh(1))/dix(2)+pom2sq(1,2)*ddu(1)+pom2sq(2,2)*ddu(2)
       end do
     case (2)
       do jj=1,4
         ddu(1)=cof(1,3)*ddx(inii(2,1),inii(jj,2),inii(2,3))+cof(2,3)*ddx(inii(2,1),inii(jj,2),inii(3,3))
         ddu(2)=cof(1,3)*ddx(inii(3,1),inii(jj,2),inii(2,3))+cof(2,3)*ddx(inii(3,1),inii(jj,2),inii(3,3))
         hrh(1)=cof(1,3)*dvl(inii(2,1),inii(jj,2),inii(2,3))+cof(2,3)*dvl(inii(2,1),inii(jj,2),inii(3,3))+&
&         pomsq(1,3)*ddz(inii(2,1),inii(jj,2),inii(2,3))+pomsq(2,3)*ddz(inii(2,1),inii(jj,2),inii(3,3))
         hrh(2)=cof(1,3)*dvl(inii(3,1),inii(jj,2),inii(2,3))+cof(2,3)*dvl(inii(3,1),inii(jj,2),inii(3,3))+&
&         pomsq(1,3)*ddz(inii(3,1),inii(jj,2),inii(2,3))+pomsq(2,3)*ddz(inii(3,1),inii(jj,2),inii(3,3))
         hh(jj,2)=(hrh(2)-hrh(1))/dix(1)+pom2sq(1,1)*ddu(1)+pom2sq(2,1)*ddu(2)

         ddu(1)=cof(1,1)*ddz(inii(2,1),inii(jj,2),inii(2,3))+cof(2,1)*ddz(inii(3,1),inii(jj,2),inii(2,3))
         ddu(2)=cof(1,1)*ddz(inii(2,1),inii(jj,2),inii(3,3))+cof(2,1)*ddz(inii(3,1),inii(jj,2),inii(3,3))
         hrh(1)=cof(1,1)*dvl(inii(2,1),inii(jj,2),inii(2,3))+cof(2,1)*dvl(inii(3,1),inii(jj,2),inii(2,3))+&
&         pomsq(1,1)*ddx(inii(2,1),inii(jj,2),inii(2,3))+pomsq(2,1)*ddx(inii(3,1),inii(jj,2),inii(2,3))
         hrh(2)=cof(1,1)*dvl(inii(2,1),inii(jj,2),inii(3,3))+cof(2,1)*dvl(inii(3,1),inii(jj,2),inii(3,3))+&
&         pomsq(1,1)*ddx(inii(2,1),inii(jj,2),inii(3,3))+pomsq(2,1)*ddx(inii(3,1),inii(jj,2),inii(3,3))
         hh(jj,1)=(hrh(2)-hrh(1))/dix(3)+pom2sq(1,3)*ddu(1)+pom2sq(2,3)*ddu(2)
       end do
     case (3)
       do jj=1,4
         ddu(1)=cof(1,1)*ddy(inii(2,1),inii(2,2),inii(jj,3))+cof(2,1)*ddy(inii(3,1),inii(2,2),inii(jj,3))
         ddu(2)=cof(1,1)*ddy(inii(2,1),inii(3,2),inii(jj,3))+cof(2,1)*ddy(inii(3,1),inii(3,2),inii(jj,3))
         hrh(1)=cof(1,1)*dvl(inii(2,1),inii(2,2),inii(jj,3))+cof(2,1)*dvl(inii(3,1),inii(2,2),inii(jj,3))+&
&         pomsq(1,1)*ddx(inii(2,1),inii(2,2),inii(jj,3))+pomsq(2,1)*ddx(inii(3,1),inii(2,2),inii(jj,3))
         hrh(2)=cof(1,1)*dvl(inii(2,1),inii(3,2),inii(jj,3))+cof(2,1)*dvl(inii(3,1),inii(3,2),inii(jj,3))+&
&         pomsq(1,1)*ddx(inii(2,1),inii(3,2),inii(jj,3))+pomsq(2,1)*ddx(inii(3,1),inii(3,2),inii(jj,3))
         hh(jj,2)=(hrh(2)-hrh(1))/dix(2)+pom2sq(1,2)*ddu(1)+pom2sq(2,2)*ddu(2)

         ddu(1)=cof(1,2)*ddx(inii(2,1),inii(2,2),inii(jj,3))+cof(2,2)*ddx(inii(2,1),inii(3,2),inii(jj,3))
         ddu(2)=cof(1,2)*ddx(inii(3,1),inii(2,2),inii(jj,3))+cof(2,2)*ddx(inii(3,1),inii(3,2),inii(jj,3))
         hrh(1)=cof(1,2)*dvl(inii(2,1),inii(2,2),inii(jj,3))+cof(2,2)*dvl(inii(2,1),inii(3,2),inii(jj,3))+&
&         pomsq(1,2)*ddy(inii(2,1),inii(2,2),inii(jj,3))+pomsq(2,2)*ddy(inii(2,1),inii(3,2),inii(jj,3))
         hrh(2)=cof(1,2)*dvl(inii(3,1),inii(2,2),inii(jj,3))+cof(2,2)*dvl(inii(3,1),inii(3,2),inii(jj,3))+&
&         pomsq(1,2)*ddy(inii(3,1),inii(2,2),inii(jj,3))+pomsq(2,2)*ddy(inii(3,1),inii(3,2),inii(jj,3))
         hh(jj,1)=(hrh(2)-hrh(1))/dix(1)+pom2sq(1,1)*ddu(1)+pom2sq(2,1)*ddu(2)
       end do
     end select
     do jj=-2,1
       grd(jj+3)=(indx(ii)+jj)*dix(ii)
     end do

!    write(std_out,'("hh: ",/,4F16.8,/,4F16.8)') ((hh(kk,jj),kk=1,4),jj=1,2)
!    write(std_out,'("grad: ",3F16.8)') (grho(kk),kk=1,3)
!    write(std_out,'("dix: ",3F16.8)') (dix(kk),kk=1,3)
!    write(std_out,'("grd: ",4F16.8)') (grd(kk),kk=1,4)
!    write(std_out,'("inii: ",4I4)') (inii(kk,ii),kk=1,4)

     do jj=1,2

!      polynomial interpolation

       do kk=1,3
         do ll=4,kk+1,-1
           hh(ll,jj)=(hh(ll,jj)-hh(ll-1,jj))/(grd(ll)-grd(ll-1))
         end do
       end do
       lder(4)=hh(4,jj)
       do kk=3,1,-1
         lder(kk)=hh(kk,jj)+(xx(ii)-grd(kk))*lder(kk+1)
       end do
       do kk=1,2
         do ll=3,kk+1,-1
           lder(ll)=lder(ll)+(xx(ii)-grd(ll-kk))*lder(ll+1)
         end do
       end do
       nn=ii+jj
       if (nn > 3) nn=nn-3
       hrho(ii,nn)=hrho(ii,nn)+lder(2)
       hrho(nn,ii)=hrho(nn,ii)+lder(2)
     end do
   end do

!  averaging of the mixed derivations obtained in different order

   do ii=1,3
     do jj=1,3
       if (ii /= jj) hrho(ii,jj)=hrho(ii,jj)/2._dp
     end do
   end do


!  write(std_out,'("xx:",3F16.8)') (xx(ii),ii=1,3)
!  write(std_out,'("hrho: ",/,3F16.8,/,3F16.8,/,3F16.8)') &
!  & ((hrho(ii,jj),ii=1,3),jj=1,3)
!  stop
!  write(std_out,'("xx:",3F16.8)') (xx(ii),ii=1,3)
!  write(std_out,'(":GRAD pred tr ",3F16.8)') grho
!  write(std_out,'(":HESSIAN pred tr",/,3F16.8,/,3F16.8,/,3F16.8)') ((hrho(ii,jj),jj=1,3),ii=1,3)


!  Transformation back to Cart. coordonnes

   call bschg1(grho,2)
   call bschg2(hrho,2)

!  write(std_out,'("hrho: ",/,3F16.8,/,3F16.8,/,3F16.8)') &
!  & ((hrho(ii,jj),ii=1,3),jj=1,3)
!  stop

   nullify(ptddx,ptddy,ptddz,ptrho)

   if (selct==1) return

 end if

!write(51,'(":GRADv ",3F16.8)') grho
!write(52,'(":LAPv ",F16.8)') hrho(1,1)+hrho(2,2)+hrho(3,3)
!write(52,'(":HESNv ",/,3F16.8,/,3F16.8,/,3F16.8)') ((hrho(ii,jj),jj=1,3),ii=1,3)

!INTERPOLATION OF THE CORE DENSITY

 if (selct/=1) then

   if (selct==2) then
     grho(:)=0._dp
     hrho(:,:)=0._dp
     rho=0._dp
   end if

!  SEARCH OF THE NEIGHBOUR ATOMS

   if (selct /= -2) then
     iat=0
     ipos=0
   end if
   rdmin=20._dp

   do jj=1,natom
     nn=typat(jj)
     rrad_nn=rrad(corlim(nn),nn)
     rrad2_nn=rrad_nn*rrad_nn
     vw1=vv(1)-xatm(1,jj)
     vw2=vv(2)-xatm(2,jj)
     vw3=vv(3)-xatm(3,jj)

     do kk=1,nnpos

       vt1=vw1-atp(1,kk)
       vt2=vw2-atp(2,kk)
       vt3=vw3-atp(3,kk)
       rr2=vt1*vt1+vt2*vt2+vt3*vt3

!      rr=vnorm(vt,0)

!      Only contribution > rhocormin (adhoc.f90) are considered

       if (rr2 < rrad2_nn .and.(.not.((selct==-2).and.(iat==jj).and.(ipos==kk)))) then
!        if (rr /= 0.0_dp) then    ! XG020629 : never test a real number against zero (not portable)
         if (rr2 > 1.0d-28) then         ! SEARCH INDEX

           rr=sqrt(rr2)
           rr_inv=1.0_dp/rr

           if (rr < rrad(1,nn)) then
             inx=-1
           elseif (rr > rrad(ndat(nn),nn)) then
             inx=ndat(nn)
           else
!            Find the index of the radius by bissection
             inmin=1
             inmax=ndat(nn)
             inx=1
             do
               if(inmax-inmin==1)exit
               inx=(inmin+inmax)/2
               if(rr>=rrad(inx,nn))then
                 inmin=inx
               else
                 inmax=inx
               end if
             end do
             inx=inmin

!            XG020629 : old coding, slower
!            inx=0
!            do while (rr >= rrad(inx+1,nn))
!            inx=inx+1
!            end do

           end if

!          Transformation matrix radial -> cart. coord
           ss=sqrt(vt1*vt1+vt2*vt2)
!          if (ss /=0._dp) then    ! XG020629 : never test a real number against zero (not portable)
           if (ss*ss > 1.0d-28) then  ! ss non-zero
!            XG020629 : very strange : only trsf(:,1) is needed in what follows ? !
!            ss_inv=1.0_dp/ss
             trsf(1,1)=vt1*rr_inv
!            trsf(1,2)=-vt2*ss_inv
!            trsf(1,3)=vt3*vt1*rr_inv*ss_inv
             trsf(2,1)=vt2*rr_inv
!            trsf(2,2)=vt1*ss_inv
!            trsf(2,3)=vt3*vt2*rr_inv*ss_inv
             trsf(3,1)=vt3*rr_inv
!            trsf(3,2)=0._dp
!            trsf(3,3)=-ss*rr_inv
!            XG020629 Not needed
!            do  ii=1,3
!            do ll=1,3
!            ches(ii,ll)=0._dp
!            end do
!            cgrad(ii)=0._dp
!            end do
           else                      ! ss zero
             do ii=1,3
               do ll=1,3
                 trsf(ii,ll)=0._dp
               end do
               trsf(ii,4-ii)=1._dp
             end do
           end if ! ss zero or non-zero

           if (inx == -1) then   ! LEFT EXTRAPOLATION y=a*x^2+b (a<0)
             val=sp2(1,nn)*0.5_dp*rr*rr/rrad(1,nn)+crho(1,nn)-sp2(1,nn)*rrad(1,nn)*0.5_dp
             cgrad(1)=sp2(1,nn)*rr/rrad(1,nn)
             ches(1,1)=sp2(1,nn)/rrad(1,nn)
           elseif (inx == ndat(nn) ) then  ! RIGHT EXTRAPOLATION y=a*exp(b*x)
             val=rrad(ndat(nn),nn)*exp(sp2(ndat(nn),nn)*(rr-rrad(ndat(nn),nn))/crho(ndat(nn),nn))
             cgrad(1)=val*sp2(ndat(nn),nn)/crho(ndat(nn),nn)
             ches(1,1)=cgrad(1)*sp2(ndat(nn),nn)/crho(ndat(nn),nn)
           else                    ! INTERPOLATION
             uu=rrad(inx+1,nn)-rrad(inx,nn)
             uu_inv=1.0_dp/uu
             aa=(rrad(inx+1,nn)-rr)*uu_inv
             bb=(rr-rrad(inx,nn))*uu_inv
             cc=(aa*aa*aa-aa)*uu*uu*0.16666666666666666_dp
             dd=(bb*bb*bb-bb)*uu*uu*0.16666666666666666_dp
             val=aa*crho(inx,nn)+bb*crho(inx+1,nn)+cc*sp3(inx,nn)+dd*sp3(inx+1,nn)
             cgrad(1)=(crho(inx+1,nn)-crho(inx,nn))*uu_inv&
&             -(3._dp*aa*aa-1._dp)*uu*0.16666666666666666_dp*sp3(inx,nn)+&
&             (3._dp*bb*bb-1._dp)*uu*0.16666666666666666_dp*sp3(inx+1,nn)
             ches(1,1)=aa*sp3(inx,nn)+bb*sp3(inx+1,nn)

           end if     ! TRANSFORMATION TO CARTEZ. COORD.

           cgrad1_rr_inv=cgrad(1)*rr_inv
           coeff=(ches(1,1)-cgrad1_rr_inv)*rr_inv*rr_inv
           cgrad(3)=trsf(3,1)*cgrad(1)
           cgrad(2)=trsf(2,1)*cgrad(1)
           cgrad(1)=trsf(1,1)*cgrad(1)
           ches(1,1)=coeff*vt1*vt1+cgrad1_rr_inv
           ches(2,2)=coeff*vt2*vt2+cgrad1_rr_inv
           ches(3,3)=coeff*vt3*vt3+cgrad1_rr_inv
           ches(1,2)=coeff*vt1*vt2 ; ches(2,1)=coeff*vt1*vt2
           ches(1,3)=coeff*vt1*vt3 ; ches(3,1)=coeff*vt1*vt3
           ches(2,3)=coeff*vt2*vt3 ; ches(3,2)=coeff*vt2*vt3

         else                                            ! case rr==0

           val=crho(1,nn)-sp2(1,nn)*rrad(1,nn)/2._dp
           do ii=1,3
             do ll=1,3
               ches(ii,ll)=0._dp
             end do
             cgrad(ii)=0._dp
             ches(ii,ii)=sp2(1,nn)/rrad(1,nn)
           end do

         end if ! rr>0 or rr==0

         do ii=1,3
           do ll=1,3
             hrho(ii,ll)=hrho(ii,ll)+ches(ii,ll)
           end do
           grho(ii)=grho(ii)+cgrad(ii)
         end do
         rho=rho+val

       end if ! rr2< rrad_nn*rrad_nn

       if (selct==-1) then
         if (rr2 < rminl(jj)*rminl(jj) ) then
           iat=jj
           ipos=kk
           rdmin=sqrt(rr2)
         end if
       elseif (selct==-2) then
         cycle
       else
         if (rr2 < rdmin*rdmin) then
           iat=jj
           ipos=kk
           rdmin=sqrt(rr2)
         end if
       end if

     end do
   end do

 end if

!write(51,'(":GRADt ",3F16.8)') grho
!write(52,'(":LAPt ",F16.8)') hrho(1,1)+hrho(2,2)+hrho(3,3)
!write(52,'(":HESNt ",/,3F16.8,/,3F16.8,/,3F16.8)') ((hrho(ii,jj),jj=1,3),ii=1,3)

!if(abs(cumul_cpu-cumul_cpu_old)>0.499)then
!write(std_out,'(a,f7.1)' )' vgh_rho : cumul_cpu=',cumul_cpu
!cumul_cpu_old=cumul_cpu
!end if

end subroutine vgh_rho
!!***

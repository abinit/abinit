!!****m* ABINIT/m_ewald
!! NAME
!!  m_ewald
!!
!! FUNCTION
!!  This module gathers routines to compute the Ewald energy and its derivatives
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2025 ABINIT group (DCA, XG, JJC, GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_ewald

use defs_basis
use m_abicore
use m_errors
use m_splines
use m_time
use m_xmpi

use m_gtermcutoff,    only : termcutoff
use m_special_funcs,  only : abi_derfc
use m_matrix,         only : matr3inv

implicit none

private

public :: ewald    ! Compute Ewald energy and derivatives with respect to xred
public :: ewald2   ! Derivative of the Ewald energy with respect to strain.
public :: ewald9   ! Compute ewald contribution to the dynamical matrix, at a given
            ! q wavevector, including anisotropic dielectric tensor and effective charges
public :: ewald9_2D! Compute ewald contribution to the dynamical matrix, at a given
            ! q wavevector, in the case of a 2D material with an external
            ! dielectric environment


contains
!!***

!!****f* m_ewald/ewald
!!
!! NAME
!! ewald
!!
!! FUNCTION
!! Compute Ewald energy and derivatives with respect to dimensionless
!! reduced atom coordinates xred.
!!
!! INPUTS
!! gmet(3,3)=metric tensor in reciprocal space (bohr^-2)
!! natom=number of atoms in unit cell
!! ntypat=numbe of type of atoms
!! rmet(3,3)=metric tensor in real space (bohr^2)
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume (bohr^3)
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!! zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!! eew=final ewald energy in hartrees
!! grewtn(3,natom)=grads of eew wrt xred(3,natom), hartrees.
!!
!! SOURCE

subroutine ewald(eew,gmet,grewtn,gsqcut,icutcoul,natom,ngfft,nkpt,ntypat,rcut,rmet,rprimd,typat,ucvol,vcutgeo,xred,zion)

!Arguments ------------------------------------
!scalars
integer,intent(in) :: icutcoul,natom,nkpt,ntypat
real(dp),intent(in) :: gsqcut,rcut,ucvol
real(dp),intent(out) :: eew
!arrays
integer,intent(in) :: ngfft(18),typat(natom)
real(dp),intent(in) :: gmet(3,3),rmet(3,3),rprimd(3,3),xred(3,natom),vcutgeo(3),zion(ntypat)
real(dp),intent(out) :: grewtn(3,natom)

!Local variables-------------------------------
!scalars
integer  :: ia,ib,ig1,ig2,ig3,ig23,ii,ir1,ir2,ir3,newg,newr,ng,nr
real(dp) :: arg,c1i,ch,chsq,derfc_arg,direct,drdta1,drdta2,drdta3,eta,fac
real(dp) :: fraca1,fraca2,fraca3,fracb1,fracb2,fracb3,gsq,gsum,phi,phr,r1
real(dp) :: minexparg
real(dp) :: r1a1d,r2,r2a2d,r3,r3a3d,recip,reta,rmagn,rsq,sumg,summi,summr,sumr
real(dp) :: t1,term ,zcut !, gcart_para, gcart_perp
!character(len=500) :: msg
!arrays
real(dp),allocatable :: gcutoff(:)

! *************************************************************************

!This is the minimum argument of an exponential, with some safety
minexparg=log(tiny(0._dp))+five

!Add up total charge and sum of $charge^2$ in cell

chsq=0._dp
ch=0._dp
do ia=1,natom
ch=ch+zion(typat(ia))
chsq=chsq+zion(typat(ia))**2
end do

!Compute eta, the Ewald summation convergence parameter,
!for approximately optimized summations:
direct=rmet(1,1)+rmet(1,2)+rmet(1,3)+rmet(2,1)+&
& rmet(2,2)+rmet(2,3)+rmet(3,1)+rmet(3,2)+rmet(3,3)
recip=gmet(1,1)+gmet(1,2)+gmet(1,3)+gmet(2,1)+&
& gmet(2,2)+gmet(2,3)+gmet(3,1)+gmet(3,2)+gmet(3,3)
!A bias is introduced, because G-space summation scales
!better than r space summation ! Note : debugging is the most
!easier at fixed eta.
zcut=SQRT(DOT_PRODUCT(rprimd(:,3),rprimd(:,3)))/2.0_dp
if(icutcoul.eq.1) then
eta=SQRT(16.0_dp/SQRT(DOT_PRODUCT(rprimd(:,1),rprimd(:,1))))
! else if (icutcoul.eq.2) then
!   zcut=SQRT(DOT_PRODUCT(rprimd(:,3),rprimd(:,3)))/2.0_dp
!   eta=217.6_dp/zcut**2.0_dp
!   eta=1.0_dp/zcut**2.0_dp
!   eta=SQRT(16.0_dp/SQRT(DOT_PRODUCT(rprimd(:,1),rprimd(:,1))))
! else if (icutcoul.eq.2) then
!   zcut=SQRT(DOT_PRODUCT(rprimd(:,3),rprimd(:,3)))/2.0_dp
!   eta=SQRT(8.0_dp/zcut)
!   eta=SQRT(16.0_dp/SQRT(DOT_PRODUCT(rprimd(:,1),rprimd(:,1))))
! else if (icutcoul.eq.2) then
!   zcut=SQRT(DOT_PRODUCT(rprimd(:,3),rprimd(:,3)))/2.0_dp
!   eta=SQRT(8.0_dp/zcut)
else
eta=pi*200.0_dp/33.0_dp*sqrt(1.69_dp*recip/direct)
end if

!Conduct reciprocal space summations
fac=pi**2/eta
gsum=0._dp
grewtn(:,:)=0.0_dp

!Initialize Gcut-off array from m_gtermcutoff
!ABI_MALLOC(gcutoff,(ngfft(1)*ngfft(2)*ngfft(3)))
call termcutoff(gcutoff,gsqcut,icutcoul,ngfft,nkpt,rcut,rprimd,vcutgeo)

!if (icutcoul.eq.3) then
!Sum over G space, done shell after shell until all
!contributions are too small.
ng=0
do
ng=ng+1
newg=0
!   Instead of this warning that most normal users do not understand (because they are doing GS calculations, and not RF calculations),
!   one should optimize this routine. But usually this is a very small fraction of any ABINIT run.
!   if (ng > 20 .and. mod(ng,10)==0) then
!      write (msg,'(3a,I10)') "Very large box of G neighbors in ewald: you probably do not want to do this.", ch10,&
!&       " If you have a metal consider setting dipdip 0.  ng = ", ng
!      ABI_WARNING(msg)
!   end if
ii=1
do ig3=-ng,ng
do ig2=-ng,ng
do ig1=-ng,ng
!        Exclude shells previously summed over
 if(abs(ig1)==ng .or. abs(ig2)==ng .or. abs(ig3)==ng .or. ng==1 ) then

!          gsq is G dot G = |G|^2
   gsq=gmet(1,1)*dble(ig1*ig1)+gmet(2,2)*dble(ig2*ig2)+&
&           gmet(3,3)*dble(ig3*ig3)+2._dp*(gmet(2,1)*dble(ig1*ig2)+&
&           gmet(3,1)*dble(ig1*ig3)+gmet(3,2)*dble(ig3*ig2))

!          Skip g=0:
   if (gsq>1.0d-20) then
     arg=fac*gsq

!            Larger arg gives 0 contribution because of exp(-arg)
     if (arg <= -minexparg ) then
!              When any term contributes then include next shell
       newg=1

       if((abs(ig1).lt.ngfft(1)).and.&
         &(abs(ig2).lt.ngfft(2)).and.&
         &(abs(ig3).lt.ngfft(3))) then
          ig23=ngfft(1)*(abs(ig2)+ngfft(2)*(abs(ig3)))
          ii=abs(ig1)+ig23+1
          !term= ( exp(-arg) + gcutoff(ii) - 1.0_dp )/gsq
          !term=exp(-arg)/gsq*gcutoff(ii)
          !term= ( exp(-arg) + gcutoff(ii) - 1.0_dp)/gsq
          term=exp(-arg)/gsq*gcutoff(ii)
       else if (icutcoul.ne.3) then
          term=zero !exp(-arg)/gsq
       else
          term=exp(-arg)/gsq
       endif

       summr = 0.0_dp
       summi = 0.0_dp


!              XG 20180531  : the two do-loops on ia should be merged, in order to spare
!              the waste of computing twice the sin and cos.

!              Note that if reduced atomic coordinates xred drift outside
!              of unit cell (outside [0,1)) it is irrelevant in the following
!              term, which only computes a phase.
       do ia=1,natom
         arg=two_pi*(ig1*xred(1,ia)+ig2*xred(2,ia)+ig3*xred(3,ia))
!                Sum real and imaginary parts (avoid complex variables)
         summr=summr+zion(typat(ia))*cos(arg)
         summi=summi+zion(typat(ia))*sin(arg)
       end do

!              The following two checks avoid an annoying underflow error msg
       if (abs(summr)<1.d-16) summr=0.0_dp
       if (abs(summi)<1.d-16) summi=0.0_dp

!              The product of term and summr**2 or summi**2 below
!              can underflow if not for checks above
       t1=term*(summr*summr+summi*summi)
       gsum=gsum+t1

       do ia=1,natom
!                Again only phase is computed so xred may fall outside [0,1).
         arg=two_pi*(ig1*xred(1,ia)+ig2*xred(2,ia)+ig3*xred(3,ia))
         phr= cos(arg)
         phi=-sin(arg)
!                (note: do not need real part, commented out)
!                c1r=(phr*summr-phi*summi)*(term*zion(typat(ia)))
         c1i=(phi*summr+phr*summi)*(term*zion(typat(ia)))
!                compute coordinate gradients
         grewtn(1,ia)=grewtn(1,ia)-c1i*ig1
         grewtn(2,ia)=grewtn(2,ia)-c1i*ig2
         grewtn(3,ia)=grewtn(3,ia)-c1i*ig3
       end do

     end if ! End condition of not larger than -minexparg
   end if ! End skip g=0
 end if ! End triple loop over G s and associated new shell condition

end do
end do
end do

!  Check if new shell must be calculated
if (newg==0) exit

end do !  End the loop on ng (new shells). Note that there is one exit from this loop.
!endif

sumg=gsum/(two_pi*ucvol)

!Stress tensor is now computed elsewhere (ewald2) hence do not need
!length scale gradients (used to compute them here).

!normalize coordinate gradients by unit cell volume ucvol
term=-2._dp/ucvol
grewtn(:,:)=grewtn(:,:)*term
!call DSCAL(3*natom,term,grewtn,1)

!Conduct real space summations
reta=sqrt(eta)
fac=2._dp*sqrt(eta/pi)
sumr=0.0_dp

!In the following a summation is being conducted over all
!unit cells (ir1, ir2, ir3) so it is appropriate to map all
!reduced coordinates xred back into [0,1).
!
!Loop on shells in r-space as was done in g-space
nr=0
do
nr=nr+1
newr=0
!   Instead of this warning that most normal users do not understand (because they are doing GS calculations, and not RF calculations),
!   one should optimize this routine. But usually this is a very small fraction of any ABINIT run.
!   if (nr > 20 .and. mod(nr,10)==0) then
!      write (msg,'(3a,I10)') "Very large box of R neighbors in ewald: you probably do not want to do this.", ch10,&
!&       " If you have a metal consider setting dipdip 0.  nr = ", nr
!      ABI_WARNING(msg)
!   end if
!
do ir3=-nr,nr
do ir2=-nr,nr
do ir1=-nr,nr
 if( abs(ir3)==nr .or. abs(ir2)==nr .or. abs(ir1)==nr .or. nr==1 )then

   do ia=1,natom
!            Map reduced coordinate xred(mu,ia) into [0,1)
     fraca1=xred(1,ia)-aint(xred(1,ia))+0.5_dp-sign(0.5_dp,xred(1,ia))
     fraca2=xred(2,ia)-aint(xred(2,ia))+0.5_dp-sign(0.5_dp,xred(2,ia))
     fraca3=xred(3,ia)-aint(xred(3,ia))+0.5_dp-sign(0.5_dp,xred(3,ia))
     drdta1=0.0_dp
     drdta2=0.0_dp
     drdta3=0.0_dp

     do ib=1,natom
!              fraca and fracb should be precomputedi and become arrays with natom dimension.
!              Also the combination with dble(ir1), dble(ir2), dble(ir3) or fraca should be done outside of the ib loop.
       fracb1=xred(1,ib)-aint(xred(1,ib))+0.5_dp-sign(0.5_dp,xred(1,ib))
       fracb2=xred(2,ib)-aint(xred(2,ib))+0.5_dp-sign(0.5_dp,xred(2,ib))
       fracb3=xred(3,ib)-aint(xred(3,ib))+0.5_dp-sign(0.5_dp,xred(3,ib))
       r1=dble(ir1)+fracb1-fraca1
       r2=dble(ir2)+fracb2-fraca2
       r3=dble(ir3)+fracb3-fraca3
       rsq=rmet(1,1)*r1*r1+rmet(2,2)*r2*r2+rmet(3,3)*r3*r3+&
&               2.0_dp*(rmet(2,1)*r2*r1+rmet(3,2)*r3*r2+rmet(3,1)*r1*r3)

!              Avoid zero denominators in 'term':
       if (rsq>=1.0d-24) then

!                Note: erfc(8) is about 1.1e-29, so do not bother with larger arg.
!                Also: exp(-64) is about 1.6e-28, so do not bother with larger arg**2 in exp.
         term=0._dp
         if (eta*rsq<64.0_dp) then
           newr=1
           rmagn=sqrt(rsq)
           arg=reta*rmagn
!                  derfc is the real(dp) complementary error function
           derfc_arg = abi_derfc(arg)
           term=derfc_arg/rmagn
           sumr=sumr+zion(typat(ia))*zion(typat(ib))*term
           term=zion(typat(ia))*zion(typat(ib))*&
&                   (term+fac*exp(-eta*rsq))/rsq
!                  Length scale grads now handled with stress tensor in ewald2
           r1a1d=rmet(1,1)*r1+rmet(1,2)*r2+rmet(1,3)*r3
           r2a2d=rmet(2,1)*r1+rmet(2,2)*r2+rmet(2,3)*r3
           r3a3d=rmet(3,1)*r1+rmet(3,2)*r2+rmet(3,3)*r3
!                  Compute terms related to coordinate gradients
           drdta1=drdta1+term*r1a1d
           drdta2=drdta2+term*r2a2d
           drdta3=drdta3+term*r3a3d
         end if
       end if ! End avoid zero denominators in'term'
     end do ! end loop over ib:

     grewtn(1,ia)=grewtn(1,ia)+drdta1
     grewtn(2,ia)=grewtn(2,ia)+drdta2
     grewtn(3,ia)=grewtn(3,ia)+drdta3
   end do ! end loop over ia:
 end if
end do ! end triple loop over real space points and associated condition of new shell
end do
end do

!  Check if new shell must be calculated
if(newr==0) exit
end do ! End loop on nr (new shells). Note that there is an exit within the loop
!
sumr=0.5_dp*sumr
fac=pi*ch**2.0_dp/(2.0_dp*eta*ucvol)

!Finally assemble Ewald energy, eew
if(icutcoul.ne.3) then
!eew=sumg+sumr-chsq*reta/sqrt(pi)-fac
eew=sumg+sumr-chsq*reta/sqrt(pi)
else
eew=sumg+sumr-chsq*reta/sqrt(pi)-fac
end if

ABI_FREE(gcutoff)

!DEBUG
!write(std_out,*)'eew=sumg+sumr-chsq*reta/sqrt(pi)-fac'
!write(std_out,*)eew,sumg,sumr,chsq*reta/sqrt(pi),fac
!ENDDEBUG

!Length scale grads handled with stress tensor, ewald2

!Output the final values of ng and nr
! write(msg, '(a,a,i4,a,i4)' )ch10,' ewald : nr and ng are ',nr,' and ',ng
! call wrtout(std_out,msg,'COLL')

end subroutine ewald
!!***

!----------------------------------------------------------------------

!!****f* m_ewald/ewald2
!!
!! NAME
!! ewald2
!!
!! FUNCTION
!! Compute the part of the stress tensor coming from the Ewald energy
!! which is calculated by derivating the Ewald energy with respect to strain.
!! See Nielsen and Martin, Phys. Rev. B 32, 3792 (1985) [[cite:Nielsen1985a]].
!! Definition of stress tensor is $(1/ucvol)*d(Etot)/d(strain(a,b))$.
!!
!! INPUTS
!! gmet(3,3)=metric tensor in reciprocal space (bohr^-2)
!! natom=number of atoms in umit cell
!! ntypat=number of type of atoms
!! rmet(3,3)=metric tensor in real space (bohr^2) (inverse transpose of gmet)
!! rprimd(3,3)=dimensional primitive translations in real space (bohr)
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume (bohr^3)
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!! zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!! $stress(6)=(1/ucvol)*gradient$ of Ewald energy with respect to strain,
!!      in hartrees/bohr^3
!! Cartesian components of stress are provided for this symmetric
!! tensor in the order 11 22 33 32 31 21.
!!
!! SOURCE

subroutine ewald2(gmet,natom,ntypat,rmet,rprimd,stress,typat,ucvol,xred,zion)

!Arguments ------------------------------------
!scalars
integer,intent(in) :: natom,ntypat
real(dp),intent(in) :: ucvol
!arrays
integer,intent(in) :: typat(natom)
real(dp),intent(in) :: gmet(3,3),rmet(3,3),rprimd(3,3),xred(3,natom)
real(dp),intent(in) :: zion(ntypat)
real(dp),intent(out) :: stress(6)

!Local variables-------------------------------
!scalars
integer :: ia,ib,ig1,ig2,ig3,ir1,ir2,ir3,newg,newr,ng,nr
real(dp) :: arg1,arg2,arg3,ch,dderfc,derfc_arg,direct,eta,fac,fraca1
real(dp) :: fraca2,fraca3,fracb1,fracb2,fracb3,g1,g2,g3,gsq,r1,r1c,r2,r2c
real(dp) :: minexparg
real(dp) :: r3,r3c,recip,reta,rmagn,rsq,summi,summr,t1,t2,t3,t4,t5,t6,term1
real(dp) :: term2,term3,term4
!arrays
real(dp) :: gprimd(3,3),strg(6),strr(6)

! *************************************************************************

!Define dimensional reciprocal space primitive translations gprimd
!(inverse transpose of rprimd)
call matr3inv(rprimd,gprimd)

!This is the minimum argument of an exponential, with some safety
minexparg=log(tiny(0._dp))+five

!Add up total charge and sum of charge^2 in cell
ch=0._dp
do ia=1,natom
ch=ch+zion(typat(ia))
end do

!Compute eta, the Ewald summation convergence parameter,
!for approximately optimized summations:
direct=rmet(1,1)+rmet(1,2)+rmet(1,3)+rmet(2,1)+&
& rmet(2,2)+rmet(2,3)+rmet(3,1)+rmet(3,2)+rmet(3,3)
recip=gmet(1,1)+gmet(1,2)+gmet(1,3)+gmet(2,1)+&
& gmet(2,2)+gmet(2,3)+gmet(3,1)+gmet(3,2)+gmet(3,3)
!Here, a bias is introduced, because G-space summation scales
!better than r space summation !
eta=pi*200.0_dp/33.0_dp*sqrt(1.69_dp*recip/direct)

fac=pi**2/eta

!Conduct reciprocal space summations
strg(1:6)=0.0_dp

!Sum over G space, done shell after shell until all
!contributions are too small
ng=0
do
ng=ng+1
newg=0

do ig3=-ng,ng
do ig2=-ng,ng
do ig1=-ng,ng

!        Exclude shells previously summed over
 if(abs(ig1)==ng .or. abs(ig2)==ng .or. abs(ig3)==ng .or. ng==1 ) then

!          Compute Cartesian components of each G
! TODO : make this a blas call, and batch things up
   g1=gprimd(1,1)*ig1+gprimd(1,2)*ig2+gprimd(1,3)*ig3
   g2=gprimd(2,1)*ig1+gprimd(2,2)*ig2+gprimd(2,3)*ig3
   g3=gprimd(3,1)*ig1+gprimd(3,2)*ig2+gprimd(3,3)*ig3
!          Compute |G|^2 (no pi factors)
   gsq=(g1**2+g2**2+g3**2)

!          skip g=0:
   if (gsq>1.0d-20) then
     arg1=fac*gsq

!            larger arg1 gives 0 contribution because of exp(-arg1)
     if (arg1<= -minexparg) then
!              When any term contributes then include next shell
       newg=1
       term1=exp(-arg1)/arg1
       summr = 0.0_dp
       summi = 0.0_dp
       do ia=1,natom
         arg2=two_pi*(ig1*xred(1,ia)+ig2*xred(2,ia)+ig3*xred(3,ia))
!                Sum real and imaginary parts (avoid complex variables)
         summr=summr+zion(typat(ia))*cos(arg2)
         summi=summi+zion(typat(ia))*sin(arg2)
       end do

!              Avoid underflow error messages
       if (abs(summr)<1.d-16) summr=0.0_dp
       if (abs(summi)<1.d-16) summi=0.0_dp

       term2=(2._dp/gsq)*(1._dp+arg1)
       t1=term2*g1*g1-1._dp
       t2=term2*g2*g2-1._dp
       t3=term2*g3*g3-1._dp
       t4=term2*g2*g3
       t5=term2*g1*g3
       t6=term2*g1*g2
       term3=term1*(summr*summr+summi*summi)
       strg(1)=strg(1)+t1*term3
       strg(2)=strg(2)+t2*term3
       strg(3)=strg(3)+t3*term3
       strg(4)=strg(4)+t4*term3
       strg(5)=strg(5)+t5*term3
       strg(6)=strg(6)+t6*term3

     end if ! End condition not being larger than -minexparg
   end if ! End skip g=0

 end if ! End triple loop and condition of new shell
end do
end do
end do

!  Check if new shell must be calculated
if (newg==0) exit
end do ! End loop on new shell. Note that there is an "exit" instruction within the loop


!Conduct real space summations
reta=sqrt(eta)
strr(1:6)=0.0_dp

!Loop on shells in r-space as was done in g-space
nr=0
do
nr=nr+1
newr=0

do ir3=-nr,nr
do ir2=-nr,nr
do ir1=-nr,nr
 if( abs(ir3)==nr .or. abs(ir2)==nr .or. abs(ir1)==nr .or. nr==1 )then

   do ia=1,natom
!            Convert reduced atomic coordinates to [0,1)
     fraca1=xred(1,ia)-aint(xred(1,ia))+0.5_dp-sign(0.5_dp,xred(1,ia))
     fraca2=xred(2,ia)-aint(xred(2,ia))+0.5_dp-sign(0.5_dp,xred(2,ia))
     fraca3=xred(3,ia)-aint(xred(3,ia))+0.5_dp-sign(0.5_dp,xred(3,ia))
     do ib=1,natom
       fracb1=xred(1,ib)-aint(xred(1,ib))+0.5_dp-sign(0.5_dp,xred(1,ib))
       fracb2=xred(2,ib)-aint(xred(2,ib))+0.5_dp-sign(0.5_dp,xred(2,ib))
       fracb3=xred(3,ib)-aint(xred(3,ib))+0.5_dp-sign(0.5_dp,xred(3,ib))
       r1=ir1+fracb1-fraca1
       r2=ir2+fracb2-fraca2
       r3=ir3+fracb3-fraca3
!              Convert from reduced to cartesian coordinates
       r1c=rprimd(1,1)*r1+rprimd(1,2)*r2+rprimd(1,3)*r3
       r2c=rprimd(2,1)*r1+rprimd(2,2)*r2+rprimd(2,3)*r3
       r3c=rprimd(3,1)*r1+rprimd(3,2)*r2+rprimd(3,3)*r3
!              Compute |r|^2
       rsq=r1c**2+r2c**2+r3c**2
       rmagn=sqrt(rsq)

!              Avoid zero denominators in 'term':
       if (rmagn>=1.0d-12) then

!                Note: erfc(8) is about 1.1e-29, so do not bother with larger arg.
!                Also: exp(-64) is about 1.6e-28, so do not bother with larger arg**2 in exp.
         arg3=reta*rmagn
         if (arg3<8.0_dp) then
           newr=1
!                  derfc computes the complementary error function
!                  dderfc is the derivative of the complementary error function
           dderfc=(-2/sqrt(pi))*exp(-eta*rsq)
           derfc_arg = abi_derfc(arg3)
           term3=dderfc-derfc_arg/arg3
           term4=zion(typat(ia))*zion(typat(ib))*term3
           strr(1)=strr(1)+term4*r1c*r1c/rsq
           strr(2)=strr(2)+term4*r2c*r2c/rsq
           strr(3)=strr(3)+term4*r3c*r3c/rsq
           strr(4)=strr(4)+term4*r2c*r3c/rsq
           strr(5)=strr(5)+term4*r1c*r3c/rsq
           strr(6)=strr(6)+term4*r1c*r2c/rsq
         end if ! End the condition of not being to large
       end if ! End avoid zero denominator

     end do ! End loop over ib:
   end do  ! End loop over ia:

 end if ! End triple loop overs real space points, and associated new shell condition
end do
end do
end do

!  Check if new shell must be calculated
if(newr==0) exit
end do ! End loop on new shells

!Finally assemble stress tensor coming from Ewald energy, stress
!(note division by unit cell volume in accordance with definition
!found in Nielsen and Martin, Phys. Rev. B 32, 3792 (1985) [[cite:Nielsen1985a]]

fac = pi/(2._dp*ucvol*eta)
stress(1)=(0.5_dp*reta*strr(1)+fac*(strg(1)+(ch**2)))/ucvol
stress(2)=(0.5_dp*reta*strr(2)+fac*(strg(2)+(ch**2)))/ucvol
stress(3)=(0.5_dp*reta*strr(3)+fac*(strg(3)+(ch**2)))/ucvol
stress(4)=(0.5_dp*reta*strr(4)+fac*strg(4))/ucvol
stress(5)=(0.5_dp*reta*strr(5)+fac*strg(5))/ucvol
stress(6)=(0.5_dp*reta*strr(6)+fac*strg(6))/ucvol

end subroutine ewald2
!!***

!!****f* m_ewald/ewald9
!! NAME
!! ewald9
!!
!! FUNCTION
!! Compute ewald contribution to the dynamical matrix, at a given
!! q wavevector, including anisotropic dielectric tensor and effective charges
!! See Phys. Rev. B 55, 10355 (1997) [[cite:Gonze1997a]], equations (72) to (75).
!! This has been generalized to quadrupoles.
!! Delivers the left hand side of Eq.(72), possibly generalized.
!!
!! INPUTS
!! acell = lengths by which lattice vectors are multiplied
!! dielt(3,3)=dielectric tensor
!! gmet(3,3) = metric in reciprocal space.
!! gprim(3,3)=dimensionless primitive translations in reciprocal space
!! natom=number of atoms in unit cell
!! qphon(3)=phonon wavevector (same system of coordinates as the reciprocal lattice vectors)
!! rmet = metric in real space
!! rprim(3,3)=dimensionless primitive translations in real space
!! sumg0: if=1, the sum in reciprocal space must include g=0,
!!  if=0, this contribution must be skipped (q=0 singularity)
!! ucvol=unit cell volume in (whatever length scale units)**3
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!! zeff(3,3,natom)=effective charge on each atom, versus electric
!!  field and atomic displacement
!! qdrp_cart(3,3,3,natom)=Quadrupole tensor on each atom in cartesian cordinates
!! option= 0: use old implementation;
!!         1: reduce the smalest argument of the exponentials to be evaluated,
!!            set eta to 1 and skip real space sum, leads to a significant speedup
!! [dipquad] = if 1, atmfrc has been build without dipole-quadrupole part
!! [quadquad] = if 1, atmfrc has been build without quadrupole-quadrupole part
!!
!! OUTPUT
!! dyew(2,3,natom,3,natom)= Ewald part of the dynamical matrix,
!!  second energy derivative wrt xred(3,natom) in Hartrees
!! Set to zero if all(zeff == zero)
!!
!! NOTES
!! 1. The q=0 part should be subtracted, by another call to
!! the present routine, with q=0. The present routine correspond
!! to the quantity written A-bar in the explanatory notes.
!! If q=0 is asked, sumg0 should be put to 0. Otherwise, it should be put to 1.
!! 2. Because this routine can be used many times in the
!! evaluation of phonons in ppddb9, it has been
!! optimized carefully. There is still possibility
!! for improvement, by using bloking on G and R!
!! 3. There can be small numerical variations due to the
!! fact that the input dielectric tensor is usually
!! not perfectly symmetric ....
!!
!! SOURCE

subroutine ewald9(acell,dielt,dyew,gmet,gprim,natom,qphon,rmet,rprim,sumg0,ucvol,xred,zeff, qdrp_cart, &
          option, dipquad, quadquad)  ! optional

!Arguments -------------------------------
!scalars
integer,intent(in) :: natom,sumg0
integer,optional,intent(in) :: option, dipquad, quadquad
real(dp),intent(in) :: ucvol
!arrays
real(dp),intent(in) :: acell(3),dielt(3,3),gmet(3,3),gprim(3,3),qphon(3)
real(dp),intent(in) :: rmet(3,3),rprim(3,3),xred(3,natom),zeff(3,3,natom)
real(dp),intent(in) :: qdrp_cart(3,3,3,natom)
real(dp),intent(out) :: dyew(2,3,natom,3,natom)

!Local variables -------------------------
!scalars
integer,parameter :: mr=10000
integer :: ia,ib,ig1,ig2,ig3,ii,ll,kk,ir,ir1,ir2,ir3,jj
integer :: info,lwork,mu,newg,newr,ng,nr,nu,ng_expxq
integer :: ewald_option
integer :: dipquad_,quadquad_
logical :: do_quadrupole
logical, save :: firstcall = .TRUE.
real(dp),parameter :: fac=4.0_dp/3.0_dp/sqrt(pi)
real(dp),parameter :: fact2=2.0_dp/sqrt(pi)
real(dp),parameter :: y2max=64.0_dp, y2min=1.0d-24
real(dp) :: cddi,cddr,cqdi,cqdr,cqqi,cqqr,g3,g4
real(dp) :: arg1,arg2,arg3,arga,c123r,c123i,c23i,c23r,detdlt,inv_detdlt
real(dp) :: direct,eta,fact1,fact3,gsq,recip,reta,reta3,inv4eta
real(dp) :: minexparg,sigma_max
real(dp) :: term1,term2,term3,term4,term5,y2,yy,invy,invy2,derfc_yy
character(len=700) :: msg
!arrays
real(dp) :: c1i(2*mr+1),c1r(2*mr+1),c2i(2*mr+1),c2r(2*mr+1),c3i(2*mr+1)
real(dp) :: c3r(2*mr+1),cosqxred(natom),wdielt(3,3),eig_dielt(3),gpq(3),gpqfac(3,3),gpqgpq(3,3)
real(dp) :: invdlt(3,3),ircar(3),ircax(3),rr(3),sinqxred(natom)
real(dp) :: xredcar(3,natom),xredcax(3,natom),xredicar(3),xredicax(3),xx(3)
real(dp) :: gprimbyacell(3,3) !,tsec(2)
real(dp),allocatable :: dyddt(:,:,:,:,:), dydqt(:,:,:,:,:,:), dyqqt(:,:,:,:,:,:,:)
real(dp),allocatable :: work(:)
complex(dpc) :: exp2piqx(natom)
complex(dpc),allocatable :: expx1(:,:), expx2(:,:), expx3(:,:)

! *********************************************************************

! This routine is expensive so skip the calculation and return zeros if zeff == zero.
! Typically this happens when the DDB file does not contains zeff but dipdip = 1 is used (default).
if (all(zeff == zero).and.all(qdrp_cart == zero)) then
dyew = zero; return
end if
do_quadrupole = any(qdrp_cart /= zero)

! Keep track of total time spent.
!call timab(1749, 1, tsec)

! Initialize dipquad and quadquad options
dipquad_=0; if(present(dipquad)) dipquad_=dipquad
quadquad_=0; if(present(quadquad)) quadquad_=quadquad

! Deactivate real space sums for quadrupolar fields or for dipdip = -1
ewald_option = 0; if (present(option)) ewald_option = option
if (do_quadrupole.and.(dipquad_==1.or.quadquad_==1)) ewald_option = 1
!ewald_option = 0

!This is the minimum argument of an exponential, with some safety
minexparg=log(tiny(0._dp))+five
if (ewald_option == 1) minexparg=-20.0_dp

! initialize complex phase factors
do ia = 1, natom
arga = two_pi*( (qphon(1))*xred(1,ia)&
          +(qphon(2))*xred(2,ia)&
          +(qphon(3))*xred(3,ia) )
exp2piqx(ia) = exp(arga*j_dpc)
end do
ng_expxq = 1000
ABI_MALLOC(expx1, (-ng_expxq:ng_expxq, natom))
ABI_MALLOC(expx2, (-ng_expxq:ng_expxq, natom))
ABI_MALLOC(expx3, (-ng_expxq:ng_expxq, natom))
do ia = 1, natom
do ig1 = -ng_expxq, ng_expxq
expx1(ig1, ia) = exp(ig1*two_pi*xred(1,ia)*j_dpc)
expx2(ig1, ia) = exp(ig1*two_pi*xred(2,ia)*j_dpc)
expx3(ig1, ia) = exp(ig1*two_pi*xred(3,ia)*j_dpc)
end do
end do

gprimbyacell = gprim
gprimbyacell(:,1) = gprimbyacell(:,1) / acell(1)
gprimbyacell(:,2) = gprimbyacell(:,2) / acell(2)
gprimbyacell(:,3) = gprimbyacell(:,3) / acell(3)

!compute eta for approximately optimized summations:
direct=rmet(1,1)+rmet(1,2)+rmet(1,3)+rmet(2,1)+&
& rmet(2,2)+rmet(2,3)+rmet(3,1)+rmet(3,2)+rmet(3,3)
recip=gmet(1,1)+gmet(1,2)+gmet(1,3)+gmet(2,1)+&
& gmet(2,2)+gmet(2,3)+gmet(3,1)+gmet(3,2)+gmet(3,3)
eta=pi*100.0_dp/33.0_dp*sqrt(1.69_dp*recip/direct)

! Compute a material-dependent width for the Gaussians that hopefully
! will make the Ewald real-space summation unnecessary.
if (ewald_option == 1) then

wdielt(:,:)=dielt(:,:)

!Diagonalize dielectric matrix
lwork=-1
ABI_MALLOC(work,(10))
call dsyev('N','U',3, wdielt, 3, eig_dielt, work, lwork,info)
lwork=nint(work(1))
ABI_FREE(work)

ABI_MALLOC(work,(lwork))
call dsyev('V','U',3, wdielt, 3, eig_dielt, work, lwork,info)
ABI_FREE(work)

!This is a tentative maximum value for the gaussian width in real space
sigma_max=three

!Set eta taking into account that the eps_inf is used as a metric in
!reciprocal space
eta=sqrt(maxval(eig_dielt))/sigma_max

if (firstcall) then
firstcall = .FALSE.
write(msg, '(4a,f9.4,9a)' ) ch10,&
' Warning : due to the use of quadrupolar fields, the width of the reciprocal space gaussians', ch10, &
' in ewald9 has been set to eta= ', eta, ' 1/bohr and the real-space sums have been neglected.', ch10, &
' One should check whether this choice leads to correct results for the specific system under study', &
' and q-point grid.',ch10, &
' It is recommended to check that calculations with dipdip=1 and -1 (both with dipquad=0 and quadquad=0)', ch10, &
' lead to identical results. Otherwise increase the resolution of the q-point grid and repeat this test.', ch10
call wrtout([ab_out,std_out], msg)
end if

!Internally eta is the square of the gaussians width
eta=eta*eta
end if

inv4eta = one / four / eta

ABI_MALLOC(dyddt,(2,3,natom,3,natom))
ABI_MALLOC(dydqt,(2,3,natom,3,natom,3))
ABI_MALLOC(dyqqt,(2,3,natom,3,natom,3,3))

dyddt = zero
dydqt = zero
dyqqt = zero

!Sum terms over g space:
ng=0
do
ng=ng+1

! if needed, update the complex phases for larger G vectors
if (ng > ng_expxq) then
!write(std_out,*)"have to realloc"
ABI_FREE(expx1)
ABI_FREE(expx2)
ABI_FREE(expx3)

ng_expxq = ng_expxq*2
! TODO: half of this space is not needed, as it contains the complex conjugate of the other half.
! present duplication avoids if statements inside the loop, however
ABI_MALLOC(expx1, (-ng_expxq:ng_expxq, natom))
ABI_MALLOC(expx2, (-ng_expxq:ng_expxq, natom))
ABI_MALLOC(expx3, (-ng_expxq:ng_expxq, natom))
do ia = 1, natom
do ig1 = -ng_expxq, ng_expxq
 expx1(ig1, ia) = exp(ig1*two_pi*xred(1,ia)*j_dpc)
 expx2(ig1, ia) = exp(ig1*two_pi*xred(2,ia)*j_dpc)
 expx3(ig1, ia) = exp(ig1*two_pi*xred(3,ia)*j_dpc)
end do
end do
end if

newg=0
do ig3=-ng,ng
do ig2=-ng,ng
do ig1=-ng,ng
 if(abs(ig1)==ng .or. abs(ig2)==ng .or. abs(ig3)==ng .or. ng==1 )then

   gpq(1)=(ig1+qphon(1))*gprimbyacell(1,1)+(ig2+qphon(2))*&
&           gprimbyacell(1,2)+(ig3+qphon(3))*gprimbyacell(1,3)
   gpq(2)=(ig1+qphon(1))*gprimbyacell(2,1)+(ig2+qphon(2))*&
&           gprimbyacell(2,2)+(ig3+qphon(3))*gprimbyacell(2,3)
   gpq(3)=(ig1+qphon(1))*gprimbyacell(3,1)+(ig2+qphon(2))*&
&           gprimbyacell(3,2)+(ig3+qphon(3))*gprimbyacell(3,3)
   gsq=zero
   do jj=1,3
     do ii=1,3
       gpqgpq(ii,jj)=gpq(ii)*gpq(jj)
       gsq=gsq+gpqgpq(ii,jj)*dielt(ii,jj)
     end do
   end do

!          Skip q=0:
   if (gsq<1.0d-20) then
     if (sumg0==1) then
       write(msg,'(a,a,a,a,a)' )&
       'The phonon wavelength should not be zero :',ch10,&
       'there are non-analytical terms that cannot be treated.',ch10,&
       'Action: subtract this wavelength from the input file.'
       ABI_ERROR(msg)
     end if

   else

     arg1=(two_pi**2)*gsq* inv4eta

!            Larger arg gives 0 contribution:
     if (arg1<= -minexparg ) then
       newg=1

!              Here calculate the term
       term1=exp(-arg1)/gsq
       do jj=1,3
         do ii=1,3
           gpqfac(ii,jj)=gpqgpq(ii,jj)*term1
         end do
       end do

! MJV: replaced old calls to cos and sin. Checked for 10 tests in v2 that max error is about 6.e-15, usually < 2.e-15
       do ia=1,natom
         cosqxred(ia)= real(exp2piqx(ia)*expx1(ig1, ia)*expx2(ig2, ia)*expx3(ig3, ia))
         sinqxred(ia)=aimag(exp2piqx(ia)*expx1(ig1, ia)*expx2(ig2, ia)*expx3(ig3, ia))
       end do

!              First, the diagonal terms
       do nu=1,3
         do ia=1,natom
           do mu=nu,3
             dyddt(1,mu,ia,nu,ia)=dyddt(1,mu,ia,nu,ia)+gpqfac(mu,nu)
           end do
         end do
       end do

!              Then, the non-diagonal ones
       do ib=2,natom
         do ia=1,ib-1
           ! phase factor dipole-dipole
           cddr=cosqxred(ia)*cosqxred(ib)+sinqxred(ia)*sinqxred(ib)
           cddi=sinqxred(ia)*cosqxred(ib)-cosqxred(ia)*sinqxred(ib)

           ! Dipole-dipole contribution
           do nu=1,3
             do mu=nu,3
               dyddt(1,mu,ia,nu,ib)=dyddt(1,mu,ia,nu,ib)+gpqfac(mu,nu)*cddr
               dyddt(2,mu,ia,nu,ib)=dyddt(2,mu,ia,nu,ib)+gpqfac(mu,nu)*cddi
             end do
           end do
         end do
       end do

       if (do_quadrupole) then
         do ib=1,natom
           do ia=1,natom

             ! phase factor for dipole-quadrupole
             cqdr=cosqxred(ia)*sinqxred(ib)-sinqxred(ia)*cosqxred(ib)
             cqdi=cosqxred(ia)*cosqxred(ib)+sinqxred(ia)*sinqxred(ib)

             ! phase factor quadrupole-quadrupole
             cqqr=cosqxred(ia)*cosqxred(ib)+sinqxred(ia)*sinqxred(ib)
             cqqi=sinqxred(ia)*cosqxred(ib)-cosqxred(ia)*sinqxred(ib)

             ! Dipole-quadrupole contribution
             do ii=1,3
               do jj=1,3
                 do kk=1,3
                   g3=gpq(ii)*gpq(jj)*gpq(kk)
                   dydqt(1,ii,ia,jj,ib,kk)=dydqt(1,ii,ia,jj,ib,kk)+g3*term1*cqdr
                   dydqt(2,ii,ia,jj,ib,kk)=dydqt(2,ii,ia,jj,ib,kk)+g3*term1*cqdi
                 end do ! kk
               end do ! jj
             end do ! ii

             ! Quadrupole-quadrupole contribution
             do ii=1,3
               do jj=1,3
                 do kk=1,3
                   do ll=1,3
                     g4 = gpq(ii)*gpq(jj)*gpq(kk)*gpq(ll)
                     dyqqt(1,ii,ia,jj,ib,kk,ll)=dyqqt(1,ii,ia,jj,ib,kk,ll)+g4*term1*cqqr
                     dyqqt(2,ii,ia,jj,ib,kk,ll)=dyqqt(2,ii,ia,jj,ib,kk,ll)+g4*term1*cqqi
                   end do
                 end do ! kk
               end do ! jj
             end do ! ii
           end do ! ia
         end do ! ib
       end if

     end if ! endif exp() argument is smaller than -minexparg
   end if ! Endif g/=0 :
 end if ! End triple summation over Gs:
end do
end do
end do

!  Check if new shell must be calculated
if(newg==0)exit
end do

!Multiplies by common factor
fact1=4.0_dp*pi/ucvol
do ib=1,natom
do ia=1,ib
do nu=1,3
do mu=nu,3
 dyddt(1,mu,ia,nu,ib)=dyddt(1,mu,ia,nu,ib)*fact1
 dyddt(2,mu,ia,nu,ib)=dyddt(2,mu,ia,nu,ib)*fact1
end do
end do
end do
end do
if (do_quadrupole) then
dydqt=dydqt*fact1/two  * two_pi
dyqqt=dyqqt*fact1/four * two_pi ** 2
end if

reta=sqrt(eta)
reta3=-eta*reta

!Calculating the inverse (transpose) of the dielectric tensor
call matr3inv(dielt,invdlt)
!Calculating the determinant of the dielectric tensor
detdlt=dielt(1,1)*dielt(2,2)*dielt(3,3)+dielt(1,3)*dielt(2,1)*&
& dielt(3,2)+dielt(1,2)*dielt(2,3)*dielt(3,1)-dielt(1,3)*&
& dielt(2,2)*dielt(3,1)-dielt(1,1)*dielt(2,3)*dielt(3,2)-&
& dielt(1,2)*dielt(2,1)*dielt(3,3)

if(detdlt<tol6)then
write(msg, '(a,es16.6,11a)' )&
'The determinant of the dielectrix matrix, detdlt=',detdlt,' is smaller than 1.0d-6.',ch10,&
'The use of the dipole-dipole model for interatomic force constants is not possible.',ch10,&
'It is likely that you have not treated the electric field perturbations,',ch10,&
'because you not are dealing with an insulator, so that',ch10,&
'your dielectric matrix was simply set to zero in the Derivative DataBase.',ch10,&
'Action: set the input variable dipdip to 0 .'
ABI_ERROR(msg)
end if

inv_detdlt = one / sqrt(detdlt)
fact3=reta3 * inv_detdlt

if (ewald_option /= 1) then
! Preparing the loop on real space
do ia=1,natom
do ii=1,3
xredcar(ii,ia)=(xred(1,ia)*acell(1)*rprim(ii,1)+&
             xred(2,ia)*acell(2)*rprim(ii,2)+&
             xred(3,ia)*acell(3)*rprim(ii,3) )*reta
end do
end do
do ia=1,natom
do ii=1,3
xredcax(ii,ia)= invdlt(1,ii)*xredcar(ii,ia)+&
             invdlt(2,ii)*xredcar(ii,ia)+&
             invdlt(3,ii)*xredcar(ii,ia)
end do
end do

! Prepare the evaluation of exp(iq*R)
do ir=-mr,mr
arg1=-two_pi*qphon(1)*ir
arg2=-two_pi*qphon(2)*ir
arg3=-two_pi*qphon(3)*ir
c1r(ir+mr+1)=cos(arg1)
c1i(ir+mr+1)=sin(arg1)
c2r(ir+mr+1)=cos(arg2)
c2i(ir+mr+1)=sin(arg2)
c3r(ir+mr+1)=cos(arg3)
c3i(ir+mr+1)=sin(arg3)
end do

do nr=1,mr
newr=0

! Begin big loop on real space vectors
do ir3=-nr,nr
do ir2=-nr,nr

! Here, construct the cosine and sine of q*R for components 2 and 3
c23r = c2r(ir2+mr+1) * c3r(ir3+mr+1) - c2i(ir2+mr+1) * c3i(ir3+mr+1)
c23i = c2i(ir2+mr+1) * c3r(ir3+mr+1) + c2r(ir2+mr+1) * c3i(ir3+mr+1)

! Also multiplies by fact3, because it is a rather economical place to do so
c23r=c23r * fact3
c23i=c23i * fact3

do ir1=-nr,nr
 if( abs(ir3)==nr .or. abs(ir2)==nr .or. abs(ir1)==nr .or. nr==1 )then

   ! This is the real part and imaginary part of the phase factor exp(iq*R)
   c123r = c1r(ir1+mr+1) * c23r - c1i(ir1+mr+1) * c23i
   c123i = c1i(ir1+mr+1) * c23r + c1r(ir1+mr+1) * c23i

   do ii=1,3
     ircar(ii)= ( ir1*acell(1)*rprim(ii,1)+&
                  ir2*acell(2)*rprim(ii,2)+&
                  ir3*acell(3)*rprim(ii,3) ) * reta
   end do
   do ii=1,3
     ircax(ii)= invdlt(1,ii)*ircar(ii)+&
                invdlt(2,ii)*ircar(ii)+&
                invdlt(3,ii)*ircar(ii)
   end do

   ! Here loops on atoms
   do ib=1,natom
     do ii=1,3
       xredicar(ii)=ircar(ii)-xredcar(ii,ib)
       xredicax(ii)=ircax(ii)-xredcax(ii,ib)
     end do
     do ia=1,ib
       do ii=1,3
         rr(ii)=xredicar(ii)+xredcar(ii,ia)
         xx(ii)=xredicax(ii)+xredcax(ii,ia)
       end do

       y2=rr(1)*xx(1)+rr(2)*xx(2)+rr(3)*xx(3)

       ! The atoms should not be too far of each other
       if (y2 < y2max) then
       ! Note: erfc(8) is about 1.1e-29, so dont bother with larger y.
       ! Also: exp(-64) is about 1.6e-28, do dont bother with larger y**2 in exp.

         ! Avoid zero denominators in term:
         if (y2 >= y2min) then
           newr=1
           yy=sqrt(y2)
           invy=1.0_dp/yy
           invy2=invy**2
           derfc_yy = abi_derfc(yy)
           term2=derfc_yy*invy*invy2
           term3=fact2*exp(-y2)*invy2
           term4=-(term2+term3)
           term5=(3.0_dp*term2+term3*(3.0_dp+2.0_dp*y2))*invy2
           do nu=1,3
             do mu=nu,3
               dyddt(1,mu,ia,nu,ib)=dyddt(1,mu,ia,nu,ib)+c123r*(xx(nu)*xx(mu)*term5+term4*invdlt(nu,mu))
               dyddt(2,mu,ia,nu,ib)=dyddt(2,mu,ia,nu,ib)+c123i*(xx(nu)*xx(mu)*term5+term4*invdlt(nu,mu))
             end do
           end do
         else
           ! If zero denominator, the atoms should be identical
           if (ia/=ib)then
             write(msg, '(5a,i0,a,i0,a)' )&
               'The distance between two atoms seem to vanish.',ch10,&
               'This is not allowed.',ch10,&
               'Action: check the input for the atoms number',ia,' and',ib,'.'
             ABI_ERROR(msg)
           else
             ! This is the correction when the atoms are identical
             do nu=1,3
               do mu=1,3
                 dyddt(1,mu,ia,nu,ib)=dyddt(1,mu,ia,nu,ib)+&
                          fac*reta3*invdlt(nu,mu) * inv_detdlt
               end do
             end do
           end if
         end if ! End the condition for avoiding zero denominators
       end if ! End the condition of too large distance between atoms
     end do
   end do ! End loop over ia and ib :
 end if ! End triple loop over real space points:
end do ! ir1
end do ! ir2
end do ! ir3

! Check if new shell must be calculated
if(newr==0)exit
if(newr==1 .and. nr==mr) ABI_BUG('mr is too small')
end do
end if ! check if should compute real part

!Now, symmetrizes
do ib=1,natom-1
do nu=1,3
do ia=ib+1,natom
do mu=nu,3
 dyddt(1,mu,ia,nu,ib)= dyddt(1,mu,ib,nu,ia)
 dyddt(2,mu,ia,nu,ib)=-dyddt(2,mu,ib,nu,ia)
end do
end do
end do
end do

do ib=1,natom
do nu=2,3
do ia=1,natom
do mu=1,nu-1
 dyddt(1,mu,ia,nu,ib)=dyddt(1,nu,ia,mu,ib)
 dyddt(2,mu,ia,nu,ib)=dyddt(2,nu,ia,mu,ib)
end do
end do
end do
end do

!Tests
!write(std_out,*)' ewald9 : take into account the effective charges '
dyew = zero
do ib=1,natom
do nu=1,3
do ia=1,natom
do mu=1,3
 do ii=1,3
   do jj=1,3
     ! dipole-dipole correction
     dyew(1,mu,ia,nu,ib)=dyew(1,mu,ia,nu,ib) + &
      zeff(ii,mu,ia)*zeff(jj,nu,ib)*dyddt(1,ii,ia,jj,ib)
     dyew(2,mu,ia,nu,ib)=dyew(2,mu,ia,nu,ib) + &
      zeff(ii,mu,ia)*zeff(jj,nu,ib)*dyddt(2,ii,ia,jj,ib)
     if (do_quadrupole) then
       do kk=1,3
         if (dipquad_==1) then
           ! dipole-quadrupole correction
           dyew(1,mu,ia,nu,ib)=dyew(1,mu,ia,nu,ib) + &
             (zeff(ii,nu,ib)*qdrp_cart(kk,jj,mu,ia) - &
              zeff(ii,mu,ia)*qdrp_cart(kk,jj,nu,ib)) * dydqt(1,ii,ia,jj,ib,kk)
           dyew(2,mu,ia,nu,ib)=dyew(2,mu,ia,nu,ib) + &
             (zeff(ii,nu,ib)*qdrp_cart(kk,jj,mu,ia) - &
              zeff(ii,mu,ia)*qdrp_cart(kk,jj,nu,ib)) * dydqt(2,ii,ia,jj,ib,kk)
         end if

         ! quadrupole-quadrupole correction
         if (quadquad_==1) then
           do ll=1,3
             dyew(1,mu,ia,nu,ib)=dyew(1,mu,ia,nu,ib) + &
             (qdrp_cart(ll,ii,mu,ia)*qdrp_cart(kk,jj,nu,ib)) * dyqqt(1,ii,ia,jj,ib,kk,ll)
             dyew(2,mu,ia,nu,ib)=dyew(2,mu,ia,nu,ib) + &
             (qdrp_cart(ll,ii,mu,ia)*qdrp_cart(kk,jj,nu,ib)) * dyqqt(2,ii,ia,jj,ib,kk,ll)
           end do
         end if
       end do
     end if

   end do
 end do
end do
end do
end do
end do

ABI_FREE(expx1)
ABI_FREE(expx2)
ABI_FREE(expx3)
ABI_FREE(dyddt)
ABI_FREE(dydqt)
ABI_FREE(dyqqt)

!call timab(1749, 2, tsec)

end subroutine ewald9
!!***

!!****f* m_ewald/ewald9_2D
!!
!! NAME
!! ewald9_2D
!!
!! FUNCTION
!! Compute the long-range electrostatics contribution to interatomic force constants
!! in the bi-dimensional (2D) case, considering the 2D is embedded in a dielectric environment
!! and has a given dielectric thickness. The singularity of the Coulomb potential is treated
!! using the Ewald summation approach. It is possible to input a more complicate model
!! with two consecutive dielectric slabs. 
!!
!! INPUTS
!! natom=number of atoms in unit cell
!! acell(3)=length of unit cell vectors
!! xred(3,natom)=reduced coordinates of the atoms
!! rprim(3,3)=unit cell vectors (unscaled)
!! dielt(3,3)=dielectric tensor of the 2D 
!! dyew(2,3,natom,3,natom)=long-range electrostatics IFCs following Ewald
!! qphon(3)=phonon wavevector in reduced coordinates
!! zeff(3,3,natom)=Born effective charge tensor
!! qdrp_cart(3,3,3,natom)=Dynamical quadrupoles
!! dielt_env=dielectric constant of the embedding environment (1 in vacuum)
!! thick(2)=dielectric thicknesses of the slab, first value correspond to the outer
!! dielectric, second to the inner dielectric slab (if any)
!! dim_msr=dimensionality of the system (indicates axis without periodicity)
!!
!! OUTPUT
!! dyew(2,3,natom,3,natom)=long-range electrostatics IFCs following Ewald
!!
!! SOURCE

subroutine ewald9_2D(natom,acell,xred,rprim,dielt,dyew,qphon,zeff,qdrp_cart,dielt_env,thick,dim_msr)  

!Arguments -------------------------------
!scalars
real(dp), intent(in) :: dielt_env
integer :: natom, dim_msr
!arrays
real(dp),intent(in) :: acell(3),thick(2),xred(3,natom),dielt(3,3),qphon(3)
real(dp),intent(in) :: rprim(3,3),zeff(3,3,natom),qdrp_cart(3,3,3,natom)
real(dp),intent(out) :: dyew(2,3,natom,3,natom)
character(len=700) :: msg

!Local variables -------------------------
!scalars
integer :: gmax,idir1,idir2,idir3,idir4,ibz1,ibz2,ibz3,ipert1,ipert2,inner_thick,ndir,mdir, rmax, rmax2
real(dp) :: detdlt, delta_perp, lambda, dielt_perp,dielt_perp1,dielt_perp2, eta,eta1,xi, dielt_eff
real(dp) :: dielt_eff1,dielt_eff2, norm_kvec, phi, qdrp_ctrcted, qdrp_ctrcted2
real(dp) :: ewald_fun, ewald_fun1, ewald_fun2, norm_real, rflct_coeff, out_thic, out_thick
real(dp) :: rflct_coeff1, rflct_coeff2, rprimd_perp, gprimd_perp, inv_qdrp, inv_qdrp2
real(dp) :: fac_erfc, fac_ewald1, fac_ewald2, fac_ewald2b,fac_exp, fac_exp1, trans_fun, fac_gauss, fac_mirror, fac_mirror1, fac_real
real(dp) :: mean2_perp, mirror_diff, mirror_parapara, mirror_paraperp, mirror_perpperp, mirror_perpperp1
real(dp) :: rvec_norm, sqrt_norm, ucsurf, xmean, sign_dip, sign_dip2
logical, save :: firstcall = .TRUE.
!arrays
integer :: periodic_dir(3)
real(dp) :: dyew_real(2,3,natom,3,natom),dyew_rec(2,3,natom,3,natom),kvec(2),invdlt_para(2,2)
real(dp) :: rprimd_para(2,2), gprimd_para(2,2), gvec(2),kvec_dielt(2,2)
real(dp) :: norm_dielt(2), dielt_para(2,2), invdlt(3,3), inv2_qdrp(3,3), inv2_qdrp2(3,3),qvec(3), qvec_para(2)
real(dp) :: diff_xcart(3), diff_xcart2(3),xcart_para(2,natom),xcart_perp(natom),fun_real(3), fun_real2(3)
real(dp) :: kvec_para(2),zeff_para(2,3,natom), zeff_perp(3,natom)
real(dp) :: xcart(3,natom), rprimd(3,3), gprimd(3,3), rvec_dielt(3)
real(dp) :: qdrp_parapara(2,2,3,natom),qdrp_perpperp(3,natom), qdrp_paraperp(2,3,natom) 
real(dp) :: rho_gerade1(2),rho_gerade2(2), rho_ungerade1(2), rho_ungerade2(2)

inner_thick = thick(1)
out_thick = thick(2)

periodic_dir(:) = 0
if (dim_msr ==2) then ! 2D along x
        periodic_dir(2) =1 ; periodic_dir(3) = 1
elseif (dim_msr==3) then ! 2D along y
        periodic_dir(1) = 1 ; periodic_dir(3) = 1
elseif (dim_msr==4) then ! 2D along z
        periodic_dir(1) = 1 ; periodic_dir(2) = 1
end if

rprimd=zero
do idir1=1,3
do idir2=1,3
rprimd(idir1,idir2) = rprim(idir1,idir2)*acell(idir2)
end do
end do


xcart=zero
do ipert1=1,natom
xcart(:,ipert1)=matmul(rprimd, xred(:,ipert1))
end do

call matr3inv(rprimd,gprimd)
zeff_para=zero ; zeff_perp=zero
qvec(:) = qphon(1)*gprimd(:,1)+qphon(2)*gprimd(:,2)
ndir=0
do idir1=1,3
if (periodic_dir(idir1)==1) then
        ndir=ndir+1
        qvec_para(ndir)=qvec(idir1)
        zeff_para(ndir,:,:) = zeff(idir1,:,:)
        xcart_para(ndir,:) = xcart(idir1,:)
else
        zeff_perp(:,:) = zeff(idir1,:,:)   
        xcart_perp(:) = xcart(idir1,:)   
        if (qvec(idir1)>tol6) then
                write(msg, '(a,es16.6,5a)')&
                        'The phonon wavevector along the confined direction is',qvec(idir1),' 1/Bohr >1.0d-6',ch10,&
                        'The phonon wavevector should be purely along the periodic direction', ch10, &
                        'when using Ewald summation in 2D. Please check your input file and structure'
                ABI_ERROR(msg)  
        end if      
end if        
end do
xmean = sum(xcart_perp)/natom
xcart_perp(:)=xcart_perp(:)-xmean

!Calculating the inverse (transpose) of the dielectric tensor
call matr3inv(dielt,invdlt)
!Calculating the determinant of the dielectric tensor
detdlt=dielt(1,1)*dielt(2,2)*dielt(3,3)+dielt(1,3)*dielt(2,1)*&
        & dielt(3,2)+dielt(1,2)*dielt(2,3)*dielt(3,1)-dielt(1,3)*&
        & dielt(2,2)*dielt(3,1)-dielt(1,1)*dielt(2,3)*dielt(3,2)-&
        & dielt(1,2)*dielt(2,1)*dielt(3,3)

if(detdlt<tol6)then
        write(msg, '(a,es16.6,11a)' )&
                'The determinant of the dielectrix matrix, detdlt=',detdlt,' is smaller than 1.0d-6.',ch10,&
                'The use of the dipole-dipole model for interatomic force constants is not possible.',ch10,&
                'It is likely that you have not treated the electric field perturbations,',ch10,&
                'because you not are dealing with an insulator, so that',ch10,&
                'your dielectric matrix was simply set to zero in the Derivative DataBase.',ch10,&
                'Action: set the input variable dipdip to 0 .'
        ABI_ERROR(msg)
end if

! The dielectric tensor must be diagonal in the confined direction
rprimd_para = zero ; rprimd_perp = zero
dielt_para = zero ; dielt_perp = zero
invdlt_para = zero
ndir=0
do idir1=1,3
if (periodic_dir(idir1)==1) then
        ndir=ndir+1
end if
mdir=0
do idir2=1,3
if (periodic_dir(idir2)==1) then
        mdir=mdir+1
end if
if ((periodic_dir(idir1)==1 .and. periodic_dir(idir2)==0) .or. &
        (periodic_dir(idir1)==0 .and. periodic_dir(idir2)==1)) then
        if (abs(dielt(idir1,idir2))>tol6 .or. abs(rprimd(idir1,idir2))>tol6) then
                write(msg, '(7a)' )&
                        'The dielectric matrix shows off-diagonal components in the confined direction larger than 1d-6',ch10,&
                        'This is forbidden when considering the Ewald summation for 2D systems. Please check if your', ch10, &
                        'confined direction is correctly specified in the anaddb input or if the vacuum size if sufficiently', ch10, &
                        'large in the confined direction to avoid spurious interactions between unit cells'
                ABI_ERROR(msg)
        end if
        qdrp_paraperp(ndir,:,:) = qdrp_cart(idir1,idir2,:,:)
elseif (periodic_dir(idir1)==0 .and. periodic_dir(idir2)==0) then
        dielt_perp = dielt(idir1,idir2)      
        rprimd_perp = rprimd(idir1,idir2)
        qdrp_perpperp(:,:) = qdrp_cart(idir1,idir2,:,:)
else
        dielt_para(ndir,mdir) = dielt(idir1,idir2)
        rprimd_para(ndir,mdir)=rprimd(idir1,idir2)     
        invdlt_para(ndir,mdir)=invdlt(idir1,idir2) 
        qdrp_parapara(ndir,mdir,:,:) = qdrp_cart(idir1,idir2,:,:)
end if
end do
end do

! First needs to determine the Gaussian broadening intrinsic to the Ewald summation. In 2D, both the real and
! reciprocal summation are related to the complementary error function. We want to restrict the real-part to
! the first Wigner cell. We use the fact that sqrt(1-e^{-x^2}) < erf(x) < sqrt(1-e^{-4x^2/pi})
! and invert those relationships to estimate the broadening required  to restrict the real-part summation of
!the Ewald summation; here fixes the threshold to 1e-6 for contribution from later unit cells
rvec_dielt = zero
ndir=0
do idir1=1,2
norm_dielt(idir1) = dot_product(rprimd_para(idir1,:),matmul(invdlt_para(:,:),rprimd_para(idir1,:)))
end do

lambda = dsqrt(maxval(norm_dielt))/dsqrt(-two*dlog(one-(one-tol9)**2))

! Now that we have computed the value of lambda, we can compute the Ewald summation, starting from the 
! reciprocal sum. We need first to estimate the max. number of reciprocal vectors we need to consider
! for the Ewald summation. We use here a stricter tolerance than for the determination of lambda
! Note that the worst case scenario is always when considering Rka=Rk'b in this case

gmax = int(dsqrt(-two*dlog(one-(one-tol12)**2))/lambda/(two_pi/dsqrt(minval(norm_dielt))))
if (firstcall) then
  firstcall = .FALSE.
  write(msg, '(6a,f9.4,2a,i3,1a)' ) ch10,&
        ' Ewald treatment of 2D long-range electrostatics interatomic force constants', ch10,  &
        ' To restrict the real-part summation of Ewald to the first unit cell, the Gaussian broadening', ch10, &
        ' has been set to ', lambda, ' 1/Bohr. For the reciprocal sum, this corresponds to max.', ch10, &
        2*gmax-1, ' Brillouin zone repetitions in either in-plane directions' 
call wrtout([ab_out,std_out], msg)
end if

ndir=0
do idir1=1,3 
if (periodic_dir(idir1)==1) then
        ndir=ndir+1
        gprimd_para(ndir,:) = gprimd(ndir,1:2)
end if
end do
dyew_rec = zero
! If one dielectric slab model, same dielectric for both regions
! Otherwise, inner dielectric ~1 and other has been computed 
! accordingly in the anaddb driver
if (out_thick> zero) then
   dielt_perp1 = one 
else
   dielt_perp1=dielt_perp
end if
dielt_perp2 = dielt_perp
do ibz1 = -gmax,gmax
  do ibz2 = -gmax,gmax
    gvec(:) = ibz1*gprimd_para(:,1)+ibz2*gprimd_para(:,2)
    kvec(:) = gvec(:) + qvec_para(:)
    kvec(:) = kvec(:)*two_pi 
    kvec_para(:) = matmul(dielt_para,kvec)
    norm_kvec = dot_product(kvec,kvec_para)
    if (abs(norm_kvec)>tol6) then !Remove G=q=0 case
       eta = dsqrt(norm_kvec/dielt_perp)
       eta1 = dsqrt(norm_kvec/dielt_perp1)
       xi = dsqrt(norm_kvec/dielt_perp2)
       dielt_eff = dsqrt(norm_kvec*dielt_perp/dot_product(kvec,kvec))
       dielt_eff1 = dsqrt(norm_kvec*dielt_perp1/dot_product(kvec,kvec))
       dielt_eff2 = dsqrt(norm_kvec*dielt_perp2/dot_product(kvec,kvec))
       ! Reflection coefficient at the dielectric interfaces
       rflct_coeff = (dielt_eff-dielt_env)/(dielt_eff+dielt_env)
       rflct_coeff2 = (dielt_eff2-dielt_env)/(dielt_eff2+dielt_env)
       trans_fun = (one+rflct_coeff2*dexp(-xi*(inner_thick-out_thick)))
       trans_fun = trans_fun/(one-rflct_coeff2*dexp(-xi*(inner_thick-out_thick)))
       rflct_coeff1 = (dielt_eff1*trans_fun-dielt_eff2)/(dielt_eff1*trans_fun+dielt_eff2)
       ! Dipole-dipole charges prefactors
       fac_exp = rflct_coeff*dexp(-eta*inner_thick)
       fac_exp1 = rflct_coeff1*dexp(-eta1*out_thick)
       fac_mirror = two*fac_exp/(one-fac_exp**2)
       fac_mirror1 = two*fac_exp1/(one-fac_exp1**2)
       ! Ewald factor in error function
       fac_ewald1= eta*dsqrt(dielt_perp)*lambda/dsqrt(two)
       do ipert1=1,natom
         do ipert2=1,natom
           delta_perp = (xcart_perp(ipert2)-xcart_perp(ipert1))
           mean2_perp = (xcart_perp(ipert2)+xcart_perp(ipert1))
           fac_ewald2= delta_perp/lambda/sqrt(two*dielt_perp)
           fac_ewald2b= delta_perp/lambda/sqrt(two*dielt_perp1)
           ! Ewald function and derivatives (eta factorized)
           ewald_fun = half*(dexp(-eta*delta_perp)*(one-erf(fac_ewald1-fac_ewald2)))+ &
               half*(dexp(eta*delta_perp)*(one-erf(fac_ewald1+fac_ewald2)))
           ewald_fun1 = half*(-dexp(-eta*delta_perp)*(one-erf(fac_ewald1-fac_ewald2)))+ &
               half*(dexp(eta*delta_perp)*(one-erf(fac_ewald1+fac_ewald2)))
           ! For second derivative, there is in principle a Gaussian term as well
           ! However, by an appropriate choice of the electrostatic gauge (mean average
           ! potential, we can neglect it. This approximation has been validated   
           ! with respect to real-space dipoles and exact calculated points
           ewald_fun2 = half*(dexp(-eta1*delta_perp)*(one-erf(fac_ewald1-fac_ewald2b)))+ &
                half*(dexp(eta1*delta_perp)*(one-erf(fac_ewald1+fac_ewald2b)))
           ewald_fun2=-ewald_fun2*eta1**2

           ! Phase factor and mirror terms
           phi = dot_product(kvec,xcart_para(:,ipert1)-xcart_para(:,ipert2))
           mirror_parapara = fac_mirror*(dcosh(eta*mean2_perp)+fac_exp*dcosh(eta*delta_perp))
           mirror_perpperp = zero
           !if (abs(fac_mirror1)>tol6) then
           mirror_perpperp = fac_mirror*(dcosh(eta*mean2_perp)-fac_exp*dcosh(eta*delta_perp))
           !print *, fac_exp1,fac_exp1**2,dcosh(eta1*mean2_perp),dcosh(eta1*delta_perp), mirror_perpperp
           !end if
           mirror_paraperp = fac_mirror*fac_exp*dsinh(eta*delta_perp)
           mirror_diff = fac_mirror*dsinh(eta*mean2_perp)
       ! Then compute the charge prefactor
           do idir1=1,3
           do idir2=1,3
           rho_gerade1(1) = -half*dot_product(kvec,matmul(qdrp_parapara(:,:,idir1,ipert1),kvec))
           rho_gerade2(1) = -half*dot_product(kvec,matmul(qdrp_parapara(:,:,idir2,ipert2),kvec))
           rho_gerade1(1) = rho_gerade1(1)+half*eta1**2*qdrp_perpperp(idir1,ipert1)
           rho_gerade2(1) = rho_gerade2(1)+half*eta1**2*qdrp_perpperp(idir2,ipert2)
           rho_gerade1(2) = -dot_product(kvec,zeff_para(:,idir1,ipert1))
           rho_gerade2(2) = -dot_product(kvec,zeff_para(:,idir2,ipert2))
           rho_ungerade1(1) = -zeff_perp(idir1,ipert1)
           rho_ungerade1(2) = eta1*dot_product(kvec,qdrp_paraperp(:,idir1,ipert1))
           rho_ungerade2(1) = -zeff_perp(idir2,ipert2)
           rho_ungerade2(2) = eta1*dot_product(kvec,qdrp_paraperp(:,idir2,ipert2))
           ! First, add the source charges (gerade gerade)
           dyew_rec(1,idir1,ipert1,idir2,ipert2)= dyew_rec(1,idir1,ipert1,idir2,ipert2)+&
                   (rho_gerade1(1)*rho_gerade2(1)+rho_gerade1(2)*rho_gerade2(2))&
                   *ewald_fun/eta/dielt_perp*cos(phi)
           dyew_rec(1,idir1,ipert1,idir2,ipert2)= dyew_rec(1,idir1,ipert1,idir2,ipert2)-&
                   (rho_gerade1(2)*rho_gerade2(1)-rho_gerade1(1)*rho_gerade2(2))&
                   *ewald_fun/eta/dielt_perp*sin(phi)
           dyew_rec(2,idir1,ipert1,idir2,ipert2)= dyew_rec(2,idir1,ipert1,idir2,ipert2)+&
                   (rho_gerade1(1)*rho_gerade2(1)+rho_gerade1(2)*rho_gerade2(2))&
                   *ewald_fun/eta/dielt_perp*sin(phi)
           dyew_rec(2,idir1,ipert1,idir2,ipert2)= dyew_rec(2,idir1,ipert1,idir2,ipert2)+&
                   (rho_gerade1(2)*rho_gerade2(1)-rho_gerade1(1)*rho_gerade2(2))&
                   *ewald_fun/eta/dielt_perp*cos(phi)        
           ! Second, add the source charges (ungerade ungerade)
           dyew_rec(1,idir1,ipert1,idir2,ipert2)= dyew_rec(1,idir1,ipert1,idir2,ipert2)+&
                   (rho_ungerade1(1)*rho_ungerade2(1)+rho_ungerade1(2)*rho_ungerade2(2))&
                   *ewald_fun2/eta1/dielt_perp1*cos(phi)
           dyew_rec(1,idir1,ipert1,idir2,ipert2)= dyew_rec(1,idir1,ipert1,idir2,ipert2)-&
                   (rho_ungerade1(2)*rho_ungerade2(1)-rho_ungerade1(1)*rho_ungerade2(2))&
                   *ewald_fun2/eta1/dielt_perp1*sin(phi)
           dyew_rec(2,idir1,ipert1,idir2,ipert2)= dyew_rec(2,idir1,ipert1,idir2,ipert2)+&
                   (rho_ungerade1(1)*rho_ungerade2(1)+rho_ungerade1(2)*rho_ungerade2(2))&
                   *ewald_fun2/eta1/dielt_perp1*sin(phi)
           dyew_rec(2,idir1,ipert1,idir2,ipert2)= dyew_rec(2,idir1,ipert1,idir2,ipert2)+&
                   (rho_ungerade1(2)*rho_ungerade2(1)-rho_ungerade1(1)*rho_ungerade2(2))&
                   *ewald_fun2/eta1/dielt_perp1*cos(phi)
           ! Third, add the source charge (gerade ungerade)
           ! For sake of consistenty, only used when there is only one dielectric thickness
           if (out_thick > zero) then
           dyew_rec(1,idir1,ipert1,idir2,ipert2)= dyew_rec(1,idir1,ipert1,idir2,ipert2)-&
                   (rho_gerade1(1)*rho_ungerade2(1)+rho_ungerade1(1)*rho_gerade2(1)&
                   +rho_gerade1(2)*rho_ungerade2(2)+rho_ungerade1(2)*rho_gerade2(2)) &
                   *ewald_fun1/eta*cos(phi)
           dyew_rec(1,idir1,ipert1,idir2,ipert2)= dyew_rec(1,idir1,ipert1,idir2,ipert2)+&
                   (rho_gerade1(2)*rho_ungerade2(1)+rho_ungerade1(2)*rho_gerade2(1)&
                    -rho_gerade1(1)*rho_ungerade2(2)-rho_ungerade1(1)*rho_gerade2(2)) &
                   *ewald_fun1/eta*sin(phi)
           dyew_rec(2,idir1,ipert1,idir2,ipert2)= dyew_rec(2,idir1,ipert1,idir2,ipert2)-&
                   (rho_gerade1(1)*rho_ungerade2(1)+rho_ungerade1(1)*rho_gerade2(1)&
                    +rho_gerade1(2)*rho_ungerade2(2)+rho_ungerade1(2)*rho_gerade2(2)) &
                   *ewald_fun1/eta*sin(phi)
           dyew_rec(2,idir1,ipert1,idir2,ipert2)= dyew_rec(2,idir1,ipert1,idir2,ipert2)-&
                   (rho_gerade1(2)*rho_ungerade2(1)+rho_ungerade1(2)*rho_gerade2(1)&
                    -rho_gerade1(1)*rho_ungerade2(2)-rho_ungerade1(1)*rho_gerade2(2)) &
                   *ewald_fun1/eta*cos(phi)
           end if
           ! Now add the interactions with the mirror charges... para para
           dyew_rec(1,idir1,ipert1,idir2,ipert2)= dyew_rec(1,idir1,ipert1,idir2,ipert2)+&
                   (rho_gerade1(1)*rho_gerade2(1)+rho_gerade1(2)*rho_gerade2(2))&
                   *mirror_parapara/eta/dielt_perp*cos(phi)
           dyew_rec(1,idir1,ipert1,idir2,ipert2)= dyew_rec(1,idir1,ipert1,idir2,ipert2)-&
                   (rho_gerade1(2)*rho_gerade2(1)-rho_gerade1(1)*rho_gerade2(2))&
                   *mirror_parapara/eta/dielt_perp*sin(phi)
           dyew_rec(2,idir1,ipert1,idir2,ipert2)= dyew_rec(2,idir1,ipert1,idir2,ipert2)+&
                   (rho_gerade1(1)*rho_gerade2(1)+rho_gerade1(2)*rho_gerade2(2))&
                   *mirror_parapara/eta/dielt_perp*sin(phi)
           dyew_rec(2,idir1,ipert1,idir2,ipert2)= dyew_rec(2,idir1,ipert1,idir2,ipert2)+&
                   (rho_gerade1(2)*rho_gerade2(1)-rho_gerade1(1)*rho_gerade2(2))&
                   *mirror_parapara/eta/dielt_perp*cos(phi)
           ! Now with perp perp
           dyew_rec(1,idir1,ipert1,idir2,ipert2)= dyew_rec(1,idir1,ipert1,idir2,ipert2)+&
                   (rho_ungerade1(1)*rho_ungerade2(1)+rho_ungerade1(2)*rho_ungerade2(2))&
                   *mirror_perpperp*eta1/dielt_perp*cos(phi)
           dyew_rec(1,idir1,ipert1,idir2,ipert2)= dyew_rec(1,idir1,ipert1,idir2,ipert2)-&
                   (rho_ungerade1(2)*rho_ungerade2(1)-rho_ungerade1(1)*rho_ungerade2(2))&
                   *mirror_perpperp*eta1/dielt_Perp*sin(phi)
           dyew_rec(2,idir1,ipert1,idir2,ipert2)= dyew_rec(2,idir1,ipert1,idir2,ipert2)+&
                   (rho_ungerade1(1)*rho_ungerade2(1)+rho_ungerade1(2)*rho_ungerade2(2))&
                   *mirror_perpperp*eta1/dielt_perp*sin(phi)
           dyew_rec(2,idir1,ipert1,idir2,ipert2)= dyew_rec(2,idir1,ipert1,idir2,ipert2)+&
                   (rho_ungerade1(2)*rho_ungerade2(1)-rho_ungerade1(1)*rho_ungerade2(2))&
                   *mirror_perpperp*eta1/dielt_perp*cos(phi)
           ! Third, add the source charge (gerade ungerade)
           if (out_thick > zero) then
           dyew_rec(1,idir1,ipert1,idir2,ipert2)= dyew_rec(1,idir1,ipert1,idir2,ipert2)+&
                   (rho_gerade1(1)*rho_ungerade2(1)+rho_ungerade1(1)*rho_gerade2(1)&
                   +rho_gerade1(2)*rho_ungerade2(2)+rho_ungerade1(2)*rho_gerade2(2)) &
                   *mirror_paraperp/eta*cos(phi)
           dyew_rec(1,idir1,ipert1,idir2,ipert2)= dyew_rec(1,idir1,ipert1,idir2,ipert2)+&
                   (rho_gerade1(2)*rho_ungerade2(1)+rho_ungerade1(2)*rho_gerade2(1)&
                    -rho_gerade1(1)*rho_ungerade2(2)-rho_ungerade1(1)*rho_gerade2(2)) &
                   *mirror_paraperp/eta*sin(phi)
           dyew_rec(2,idir1,ipert1,idir2,ipert2)= dyew_rec(2,idir1,ipert1,idir2,ipert2)-&
                   (rho_gerade1(1)*rho_ungerade2(1)+rho_ungerade1(1)*rho_gerade2(1)&
                    +rho_gerade1(2)*rho_ungerade2(2)+rho_ungerade1(2)*rho_gerade2(2)) &
                   *mirror_paraperp/eta*sin(phi)
           dyew_rec(2,idir1,ipert1,idir2,ipert2)= dyew_rec(2,idir1,ipert1,idir2,ipert2)+&
                   (rho_gerade1(2)*rho_ungerade2(1)+rho_ungerade1(2)*rho_gerade2(1)&
                    -rho_gerade1(1)*rho_ungerade2(2)-rho_ungerade1(1)*rho_gerade2(2)) &
                   *mirror_paraperp/eta*cos(phi)
           ! Finally, there is a term on the sum of charge, only for mirror charges
                   ! Third, add the source charge (gerade ungerade)
           dyew_rec(1,idir1,ipert1,idir2,ipert2)= dyew_rec(1,idir1,ipert1,idir2,ipert2)-&
                   (rho_gerade1(1)*rho_ungerade2(1)-rho_ungerade1(1)*rho_gerade2(1)&
                   +rho_gerade1(2)*rho_ungerade2(2)-rho_ungerade1(2)*rho_gerade2(2)) &
                   *mirror_diff/eta*cos(phi)
           dyew_rec(1,idir1,ipert1,idir2,ipert2)= dyew_rec(1,idir1,ipert1,idir2,ipert2)-&
                   (rho_gerade1(2)*rho_ungerade2(1)-rho_ungerade1(2)*rho_gerade2(1)&
                    -rho_gerade1(1)*rho_ungerade2(2)+rho_ungerade1(1)*rho_gerade2(2)) &
                   *mirror_diff/eta*sin(phi)
           dyew_rec(2,idir1,ipert1,idir2,ipert2)= dyew_rec(2,idir1,ipert1,idir2,ipert2)+&
                   (rho_gerade1(1)*rho_ungerade2(1)-rho_ungerade1(1)*rho_gerade2(1)&
                    +rho_gerade1(2)*rho_ungerade2(2)-rho_ungerade1(2)*rho_gerade2(2)) &
                   *mirror_diff/eta*sin(phi)
           dyew_rec(2,idir1,ipert1,idir2,ipert2)= dyew_rec(2,idir1,ipert1,idir2,ipert2)-&
                   (rho_gerade1(2)*rho_ungerade2(1)-rho_ungerade1(2)*rho_gerade2(1)&
                    +rho_gerade1(1)*rho_ungerade2(2)+rho_ungerade1(1)*rho_gerade2(2)) &
                   *mirror_diff/eta*sin(phi)
           end if
           end do
           end do      
         end do
       end do   
     end if
  end do
end do 
ucsurf = rprimd_para(1,1)*rprimd_para(2,2)-rprimd_para(1,2)*rprimd_para(2,1)
! Renormalize by surface of periodic 2D lattice and out-of-plane dielectric constant
dyew_rec = dyew_rec*(two_pi)/ucsurf
! Reciprocal summation completes. Remains some real-space contribution from first
! unit cells (ipert1 neq ipert2), impacting all IFCs the same way
dyew_real = zero
do ipert1=1,natom
  do ipert2=1,natom
    diff_xcart(:) =xcart(:,ipert2)-xcart(:,ipert1)
    rvec_dielt(:) = matmul(invdlt,diff_xcart)
    rvec_norm = dot_product(diff_xcart,rvec_dielt)
    sqrt_norm = dsqrt(rvec_norm)
    fac_erfc = sqrt_norm/dsqrt(two)/lambda
    fac_gauss = -rvec_norm/two/lambda**2
    do idir1=1,3
      do idir2=1,3
        if (ipert1 .NE. ipert2) then ! Only off-sites contributions
          ! First, contribution from eps^-1 (Rk'b-Rka) eps^-1
          fac_real = three*(one-erf(fac_erfc))/sqrt_norm**5+6*exp(fac_gauss)/rvec_norm**2/&
          sqrt(two_pi)/lambda+two*exp(fac_gauss)/rvec_norm/sqrt(two_pi)/lambda**3
          dyew_real(1,idir1,ipert1,idir2,ipert2) = dyew_real(1,idir1,ipert1,idir2,ipert2)+&
          fac_real*rvec_dielt(idir1)*rvec_dielt(idir2) 
          ! Second contribution from esp^-1(alpha,beta)  
          fac_real = (one-erf(fac_erfc))/sqrt_norm**3+exp(fac_gauss)/rvec_norm/sqrt(two_pi)/lambda 
          dyew_real(1,idir1,ipert1,idir2,ipert2)=dyew_real(1,idir1,ipert1,idir2,ipert2)-&
          fac_real*dielt(idir1,idir2)
        end if          
      end do          
    end do
  end do
end do  
dyew_real = dyew_real / dsqrt(detdlt)
dyew = dyew_real + dyew_rec
end subroutine ewald9_2D

end module m_ewald
!!***

!!****m* ABINIT/m_ewald
!! NAME
!!  m_ewald
!!
!! FUNCTION
!!  This module gathers routines to compute the Ewald energy and its derivatives
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (DCA, XG, JJC, GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
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
 use m_symtk,          only : matr3inv

 implicit none

 private

 public :: ewald    ! Compute Ewald energy and derivatives with respect to xred
 public :: ewald2   ! Derivative of the Ewald energy with respect to strain.
 public :: ewald9   ! Compute ewald contribution to the dynamical matrix, at a given
                    ! q wavevector, including anisotropic dielectric tensor and effective charges

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
!! PARENTS
!!      m_setvtr
!!
!! CHILDREN
!!      dsyev,matr3inv,timab,wrtout
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
 real(dp) :: t1,term,zcut
 !character(len=500) :: message
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
if(icutcoul.eq.1) then
   eta=SQRT(16.0_dp/SQRT(DOT_PRODUCT(rprimd(:,1),rprimd(:,1))))
 else if (icutcoul.eq.2) then
   zcut=SQRT(DOT_PRODUCT(rprimd(:,3),rprimd(:,3)))/2.0_dp
   eta=SQRT(8.0_dp/zcut)
 else
   eta=pi*200.0_dp/33.0_dp*sqrt(1.69_dp*recip/direct)
 end if

!Conduct reciprocal space summations
 fac=pi**2/eta
 gsum=0._dp
 grewtn(:,:)=0.0_dp

 !Initialize Gcut-off array from m_gtermcutoff
 !ABI_ALLOCATE(gcutoff,(ngfft(1)*ngfft(2)*ngfft(3)))
 call termcutoff(gcutoff,gsqcut,icutcoul,ngfft,nkpt,rcut,rprimd,vcutgeo)

!Sum over G space, done shell after shell until all
!contributions are too small.
 ng=0
 do
   ng=ng+1
   newg=0
!   Instead of this warning that most normal users do not understand (because they are doing GS calculations, and not RF calculations),
!   one should optimize this routine. But usually this is a very small fraction of any ABINIT run.
!   if (ng > 20 .and. mod(ng,10)==0) then
!      write (message,'(3a,I10)') "Very large box of G neighbors in ewald: you probably do not want to do this.", ch10,&
!&       " If you have a metal consider setting dipdip 0.  ng = ", ng
!      MSG_WARNING(message)
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

!              The following two checks avoid an annoying underflow error message
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
!      write (message,'(3a,I10)') "Very large box of R neighbors in ewald: you probably do not want to do this.", ch10,&
!&       " If you have a metal consider setting dipdip 0.  nr = ", nr
!      MSG_WARNING(message)
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
   eew=sumg+sumr-chsq*reta/sqrt(pi)
 else
   eew=sumg+sumr-chsq*reta/sqrt(pi)-fac
 end if

 ABI_DEALLOCATE(gcutoff) 

!DEBUG
write(std_out,*)'eew=sumg+sumr-chsq*reta/sqrt(pi)-fac'
write(std_out,*)eew,sumg,sumr,chsq*reta/sqrt(pi),fac
!ENDDEBUG

!Length scale grads handled with stress tensor, ewald2

!Output the final values of ng and nr
! write(message, '(a,a,i4,a,i4)' )ch10,' ewald : nr and ng are ',nr,' and ',ng
! call wrtout(std_out,message,'COLL')

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
!! PARENTS
!!      m_stress
!!
!! CHILDREN
!!      dsyev,matr3inv,timab,wrtout
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
!! See Phys. Rev. B 55, 10355 (1997) [[cite:Gonze1997a]], equations (71) to (75).
!!
!! INPUTS
!! acell = lengths by which lattice vectors are multiplied
!! dielt(3,3)=dielectric tensor
!! gmet(3,3) = metric in reciprocal space.
!! gprim(3,3)=dimensionless primitive translations in reciprocal space
!! natom=number of atoms in unit cell
!! qphon(3)=phonon wavevector (same system of coordinates as the
!!  reciprocal lattice vectors)
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
!! PARENTS
!!      m_dynmat,m_effective_potential,m_ifc
!!
!! CHILDREN
!!      dsyev,matr3inv,timab,wrtout
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
 integer,parameter :: mr=10000,ny2_spline=1024*10
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
 character(len=700) :: message
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

 ! Deactivate real space sums for quadrupolar fileds or for dipdip=-1
 ewald_option = 0; if (present(option)) ewald_option = option
 if (do_quadrupole.and.(dipquad_==1.or.quadquad_==1)) ewald_option = 1

!This is the minimum argument of an exponential, with some safety
 minexparg=log(tiny(0._dp))+five
 if (ewald_option == 1) minexparg=-20.0_dp

! initialize complex phase factors
 do ia = 1, natom
   arga = two_pi*( (qphon(1))*xred(1,ia)&
&                 +(qphon(2))*xred(2,ia)&
&                 +(qphon(3))*xred(3,ia) )
   exp2piqx(ia) = exp(arga*j_dpc)
 end do
 ng_expxq = 1000
 ABI_ALLOCATE(expx1, (-ng_expxq:ng_expxq, natom))
 ABI_ALLOCATE(expx2, (-ng_expxq:ng_expxq, natom))
 ABI_ALLOCATE(expx3, (-ng_expxq:ng_expxq, natom))
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
 ! will make the Ewald real-space summation innecessary.
 if (ewald_option == 1) then

   wdielt(:,:)=dielt(:,:)

   !Diagonalize dielectric matrix
   lwork=-1
   ABI_ALLOCATE(work,(10))
   call dsyev('N','U',3, wdielt, 3, eig_dielt, work, lwork,info)
   lwork=nint(work(1))
   ABI_DEALLOCATE(work)

   ABI_ALLOCATE(work,(lwork))
   call dsyev('V','U',3, wdielt, 3, eig_dielt, work, lwork,info)
   ABI_DEALLOCATE(work)

   !This is a tentative maximum value for the gaussian width in real space
   sigma_max=three

   !Set eta taking into account that the eps_inf is used as a metric in
   !reciprocal space
   eta=sqrt(maxval(eig_dielt))/sigma_max

   if (firstcall) then
     firstcall = .FALSE.
     write(message, '(4a,f9.4,9a)' ) ch10,&
    &' Warning : due to the use of quadrupolar fields, the width of the reciprocal space gaussians', ch10, &
    &' in ewald9 has been set to eta= ', eta, ' 1/bohr and the real-space sums have been neglected.', ch10, &
    &' One should check whether this choice leads to correct results for the specific system under study', &
    &' and q-point grid.',ch10, &
    &' It is recommended to check that calculations with dipdip=1 and -1 (both with dipquad=0 and quadquad=0)', ch10, &
    &' lead to identical results. Otherwise increase the resolution of the q-point grid and repeat this test.', ch10
     call wrtout([ab_out,std_out],message,'COLL')
   end if

   !Internally eta is the square of the gaussians width
   eta=eta*eta

 end if

 inv4eta = one / four / eta

 ABI_ALLOCATE(dyddt,(2,3,natom,3,natom))
 ABI_ALLOCATE(dydqt,(2,3,natom,3,natom,3))
 ABI_ALLOCATE(dyqqt,(2,3,natom,3,natom,3,3))

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
     ABI_DEALLOCATE(expx1)
     ABI_DEALLOCATE(expx2)
     ABI_DEALLOCATE(expx3)

     ng_expxq = ng_expxq*2
! TODO: half of this space is not needed, as it contains the complex conjugate of the other half.
! present duplication avoids if statements inside the loop, however
     ABI_ALLOCATE(expx1, (-ng_expxq:ng_expxq, natom))
     ABI_ALLOCATE(expx2, (-ng_expxq:ng_expxq, natom))
     ABI_ALLOCATE(expx3, (-ng_expxq:ng_expxq, natom))
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
               write(message,'(a,a,a,a,a)' )&
&               'The phonon wavelength should not be zero :',ch10,&
&               'there are non-analytical terms that cannot be treated.',ch10,&
&               'Action: subtract this wavelength from the input file.'
               MSG_ERROR(message)
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
   write(message, '(a,es16.6,11a)' )&
&   'The determinant of the dielectrix matrix, detdlt=',detdlt,' is smaller than 1.0d-6.',ch10,&
&   'The use of the dipole-dipole model for interatomic force constants is not possible.',ch10,&
&   'It is likely that you have not treated the electric field perturbations,',ch10,&
&   'because you not are dealing with an insulator, so that',ch10,&
&   'your dielectric matrix was simply set to zero in the Derivative DataBase.',ch10,&
&   'Action: set the input variable dipdip to 0 .'
   MSG_ERROR(message)
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
                     write(message, '(5a,i0,a,i0,a)' )&
                       'The distance between two atoms seem to vanish.',ch10,&
                       'This is not allowed.',ch10,&
                       'Action: check the input for the atoms number',ia,' and',ib,'.'
                     MSG_ERROR(message)
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
   if(newr==1 .and. nr==mr) MSG_BUG('mr is too small')
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

 ABI_DEALLOCATE(expx1)
 ABI_DEALLOCATE(expx2)
 ABI_DEALLOCATE(expx3)
 ABI_DEALLOCATE(dyddt)
 ABI_DEALLOCATE(dydqt)
 ABI_DEALLOCATE(dyqqt)

 !call timab(1749, 2, tsec)

end subroutine ewald9
!!***

end module m_ewald
!!***

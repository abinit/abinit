!{\src2tex{textfont=tt}}
!!****f* ABINIT/entropyrec
!! NAME
!! entropyrec
!! 
!! FUNCTION
!! This routine computes the local part of the entropy at a point using a path integral, 
!! in the recursion method.
!! 
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  an, bn2 : coefficient given by the recursion.
!!  nrec=order of recursion
!!  trotter=trotter parameter
!!  multce=a multiplicator for computing entropy ; 2 for non-spin-polarized system
!!  debug_rec=debug variable
!!  n_pt_integ=number of points of integration for the path integral 
!!  xmax =max point of integration on the real axis

!! OUTPUT
!!  ent_out=entropy at the point
!!  ent_out1,ent_out2,ent_out3,ent_out4=debug entropy at the point
!!  
!! PARENTS
!!      vtorhorec
!!
!! CHILDREN
!!      timab,wrtout
!!
!! NOTES
!!  at this time :
!!       - multce should be not used
!!       - the routine should be integraly rewrited and use the routine recursion. 
!!       - only modified for p /= 0
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine entropyrec(an,bn2,nrec,trotter,ent_out,multce,debug_rec, &
&                     n_pt_integ,xmax,&
&                     ent_out1,ent_out2,ent_out3,ent_out4)


 use defs_basis
 use m_profiling_abi

 use m_time,             only : timab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'entropyrec'
 use interfaces_14_hidewrite
 use interfaces_18_timing
!End of the abilint section

 implicit none
  
!Arguments -------------------------------
!scalars
 integer,intent(in) :: n_pt_integ,nrec,trotter
 logical,intent(in) :: debug_rec
 real(dp), intent(in) :: multce,xmax
 real(dp),intent(out) :: ent_out,ent_out1,ent_out2,ent_out3,ent_out4
!arrays
 real(dp),intent(in) :: an(0:nrec),bn2(0:nrec)
 
!Local variables-------------------------------
!scalars
 integer, parameter :: level = 7
 integer, save :: first_en = 1
 integer :: ii,kk,n_pt_integ_path2,n_pt_integ_path3
 real(dp) :: arg,epsilo,step,twotrotter,xmin,dr_step
 complex(dpc) :: D,Dnew,Dold,N,Nnew,Nold,dz_path,ent_acc,ent_acc1,ent_acc2
 complex(dpc) :: ent_acc3,ent_acc4
 complex(dpc) :: funczero,z_path,zj
 complex(dpc) ::delta_calc
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp) :: iif,factor


! *************************************************************************


 call timab(610,1,tsec)
 
!structured debugging if debug_rec=T : print detailled result the first time we enter entropyrec
 
 if(debug_rec .and. first_en==1)then
   write(msg,'(a)')' ' 
   call wrtout(std_out,msg,'PERS')
   write(msg,'(a)')' entropyrec : enter ' 
   call wrtout(std_out,msg,'PERS')
   write(msg,'(a,i6)')'n_pt_integ ' , n_pt_integ
   call wrtout(std_out,msg,'COLL')
 end if
 
 ent_out = zero
 ent_out1 = zero
 ent_out2 = zero
 ent_out3 = zero
 ent_out4 = zero
 ent_acc = czero
 ent_acc1 = czero
 ent_acc2 = czero
 ent_acc3 = czero
 ent_acc4 = czero

!path parameters 
 twotrotter = max(two*real(trotter,dp),one)
 if(trotter==0)then
   factor = tol5
   arg =pi*three_quarters
   zj = cmplx(-one,one-sin(arg),dp)
 else
   factor = xmax/ten
   arg = pi/twotrotter
   zj = cmplx( cos(arg) , sin(arg),dp )
 end if

 epsilo = factor*sin( arg )
 xmin = factor*cos( arg )
 step = (xmax-xmin)/real(n_pt_integ,dp)

!####################################################################
![xmax + i*epsilo,xmin + i*epsilo]
 dr_step = one/real(n_pt_integ,dp)
 path1:  do ii = 0,n_pt_integ
   z_path = cmplx(xmin+real(ii,dp)*(xmax-xmin)*dr_step,epsilo,dp)
   dz_path = -cmplx((xmax-xmin)*dr_step,zero,dp)
   
   Nold = czero
   Dold = cone
   N = cone
   D = z_path - cmplx(an(0),zero,dp)
   
   do kk=1,nrec
     Nnew = (z_path - cmplx(an(kk),zero,dp))*N - cmplx(bn2(kk),zero,dp)*Nold
     Dnew = (z_path - cmplx(an(kk),zero,dp))*D - cmplx(bn2(kk),zero,dp)*Dold   

     Nold = N
     Dold = D
     N = Nnew
     D = Dnew

     if(kk/=nrec)then
       if((bn2(kk+1)<tol14))exit
     end if
   end do

!  <r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   delta_calc = func1_rec(z_path**twotrotter)*(N/D)*dz_path
   if(ii==0.or.ii==n_pt_integ)then
     ent_acc  = ent_acc  + half*delta_calc
     ent_acc1 = ent_acc1 + half*delta_calc
   else
     ent_acc  = ent_acc  + delta_calc
     ent_acc1 = ent_acc1 + delta_calc
   end if
 end do path1
 

!####################################################################
![1/2zj,0]
 if(epsilo/step>100.d0)then
   n_pt_integ_path2 = int((factor*abs(zj))/step)+1
 else
   n_pt_integ_path2 = 100
 end if

 if(trotter/=0)then
   n_pt_integ_path3 = 0
   dr_step = one/real(n_pt_integ_path2,dp)
   dz_path = -cmplx(xmin,epsilo,dp)*dr_step
   path5:  do ii = 0,n_pt_integ_path2
     z_path = cmplx(real(ii,dp)*xmin,real(ii,dp)*epsilo,dp)*dr_step
     if(abs(z_path)>tol14)then
       Nold = czero
       Dold = cone
       N = cone
       D = z_path - cmplx(an(0),zero,dp)
       do kk=1,nrec
         Nnew = (z_path - cmplx(an(kk),zero,dp))*N - cmplx(bn2(kk),zero,dp)*Nold
         Dnew = (z_path - cmplx(an(kk),zero,dp))*D - cmplx(bn2(kk),zero,dp)*Dold
         Nold = N
         Dold = D
         N = Nnew
         D = Dnew
         if(kk/=nrec)then
           if((bn2(kk+1)<tol14))exit
         end if
       end do

!      <r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
       if(abs(z_path)**twotrotter>tiny(one)) then
         funczero = func1_rec(z_path**twotrotter)
       else 
         funczero = czero
       end if
       delta_calc = funczero*N/D*dz_path
       if(ii==0.or.ii==n_pt_integ_path2)then
         ent_acc  = ent_acc  + half*delta_calc
         if(debug_rec) ent_acc3 = ent_acc3 + half*delta_calc
       else
         ent_acc  = ent_acc  + funczero*delta_calc
         if(debug_rec) ent_acc3 = ent_acc3 + funczero*delta_calc 
       end if
     end if
   end do path5

 else  ! trotter==0

   n_pt_integ_path3 = max(100,int((epsilo*half*pi)/real(step,dp))+1)
   dr_step = one/real(n_pt_integ_path3,dp)
   path6:  do ii = 0,n_pt_integ_path3
     iif=half*pi*real(ii,dp)*dr_step
     z_path = epsilo*cmplx(-cos(iif),1-sin(iif),dp)
     dz_path = epsilo*cmplx(sin(iif),-cos(iif),dp)*half*pi*dr_step
     if(abs(z_path)**twotrotter>tol14)then
       Nold = czero
       Dold = cone
       N = cone
       D = z_path - cmplx(an(0),zero,dp)
       do kk=1,nrec
         Nnew = (z_path - cmplx(an(kk),zero,dp))*N - cmplx(bn2(kk),zero,dp)*Nold
         Dnew = (z_path - cmplx(an(kk),zero,dp))*D - cmplx(bn2(kk),zero,dp)*Dold
         Nold = N
         Dold = D
         N = Nnew
         D = Dnew
         if(kk/=nrec .and. bn2(kk+1)<tol14) exit !-EXIT
       end do

!      <r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
       delta_calc = func1_rec(z_path**twotrotter) * N/D * dz_path
       if(ii==0.or.ii==n_pt_integ_path3)then
         ent_acc  = ent_acc + half*delta_calc 
         if(debug_rec) ent_acc3 = ent_acc3 + half*delta_calc 
       else
         ent_acc  = ent_acc + delta_calc    !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
         if(debug_rec) ent_acc3 = ent_acc3 + delta_calc  !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
       end if
     end if
   end do path6

 end if

 if(first_en==1 .and. debug_rec) then
   write(msg,'(a,i5,2a,i5,2a,i5,2a,es11.4,2a,es11.4,2a,es11.4)')&
&   'n_pt_path  =',n_pt_integ,ch10,&
&   'n_pt_path2 =',n_pt_integ_path2,ch10,&
&   'n_pt_path3 =',n_pt_integ_path3,ch10,&
&   'xmin       =',xmin,ch10,&
&   'xmax       =',xmax,ch10,&
&   'epsilon    =',epsilo
   call wrtout(std_out,msg,'COLL')
   first_en = 0
 end if

!####################################################################
![xmax,xmax+i*epsilo]
 dr_step = one/real(n_pt_integ_path2,dp)
 dz_path = cmplx(zero,epsilo*dr_step,dp)
 path4:  do ii = 0,n_pt_integ_path2
   z_path = cmplx(xmax,real(ii,dp)*epsilo*dr_step,dp)

   Nold = czero
   Dold = cone
   N = cone
   D = z_path - cmplx(an(0),zero,dp)

   do kk=1,nrec
     Nnew = (z_path - cmplx(an(kk),zero,dp))*N - cmplx(bn2(kk),zero,dp)*Nold
     Dnew = (z_path - cmplx(an(kk),zero,dp))*D - cmplx(bn2(kk),zero,dp)*Dold

     Nold = N
     Dold = D
     N = Nnew
     D = Dnew

     if(kk/=nrec)then
       if((bn2(kk+1)<tol14))exit
     end if
   end do

!  <r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   delta_calc = func1_rec(z_path**twotrotter)*N/D*dz_path
   if(ii==0.or.ii==n_pt_integ_path2)then

     ent_acc =  ent_acc  + half*delta_calc 
     if(debug_rec) ent_acc2 = ent_acc2 + half*delta_calc 
   else
     ent_acc  = ent_acc  + delta_calc      
     if(debug_rec) ent_acc2 = ent_acc2 + delta_calc
   end if
 end do path4


 ent_out  = multce*real(ent_acc*cmplx(zero,-piinv,dp),dp)
 if(debug_rec) then
   ent_out1 = multce*real(ent_acc1*cmplx(zero,-piinv,dp),dp)
   ent_out2 = multce*real(ent_acc2*cmplx(zero,-piinv,dp),dp)
   ent_out3 = multce*real(ent_acc3*cmplx(zero,-piinv,dp),dp)
   ent_out4 = multce*real(ent_acc4*cmplx(zero,-piinv,dp),dp)
 end if

 call timab(610,2,tsec)

 contains

!function to integrate over the path
!func1_rec(z_path,twotrotter) =  ( z_path**twotrotter/(1+z_path**twotrotter)*log(1+1/z_path**twotrotter)+&    !- f*ln(f)
!&1/(1+z_path**twotrotter)*log(1+z_path**twotrotter))       !- (1-f)*ln(1-f)

!func1_rec(z_path_pow) =   z_path_pow/(cone+z_path_pow)*log(cone+cone/z_path_pow)+&    !- f*ln(f)
!&cone/(cone+z_path_pow)*log(cone+z_path_pow)       !- (1-f)*ln(1-f)

!other expression of func for a path like ro(t)*exp(2*i*pi/(2*p)*(j+1/2))

   function func1_rec(z)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'func1_rec'
!End of the abilint section

   implicit none

   complex(dpc) :: func1_rec
   complex(dpc),intent(in) :: z

   func1_rec =   z/(cone+z)*log(cone+cone/z)+ cone/(cone+z)*log(cone+z)

 end function func1_rec
 
end subroutine entropyrec
!!***

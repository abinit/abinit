!{\src2tex{textfont=tt}}
!!****f* ABINIT/gran_potrec
!! NAME
!! gran_potrec
!! 
!! FUNCTION
!! This routine computes the local part of the grand-potential at a point using a path integral, 
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
!!  mult=a multiplicator for computing grand-potential (2 for non-spin-polarized system)
!!  debug_rec=debugging variable
!!  n_pt_integ=points for computation of path integral  
!!  xmax= maximum point on the x-axis for integration
!! 
!! OUTPUT
!!  ene_out=grand-potential at the point
!!  if debug_rec=T then ene_out1,ene_out2,ene_out3,ene_out4 are
!!  the different path branch contriubutions to the grand-potential.
!!  In reality it is not the gren potential but the
!!  grand-potential (omega=-PV) divided by -T 
!!    
!!
!! PARENTS
!!      vtorhorec
!!
!! CHILDREN
!!      timab,wrtout
!!
!! NOTES
!!  in reality it is not the gren potential but the grand-potential (omega=-PV) divided by -T 
!!  at this time :
!!       - mult should be not used
!!       - the routine should be integraly rewrited and use the routine recursion. 
!!       - only modified for p /= 0
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine gran_potrec(an,bn2,nrec,trotter,ene_out, mult, &
&                     debug_rec,n_pt_integ,xmax,&
&                     ene_out1,ene_out2,ene_out3,ene_out4)

 use defs_basis
 use m_profiling_abi

 use m_time,        only : timab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gran_potrec'
 use interfaces_14_hidewrite
 use interfaces_18_timing
!End of the abilint section

 implicit none
  
!Arguments -------------------------------
!scalars
 integer,intent(in) :: n_pt_integ,nrec,trotter
 logical,intent(in) :: debug_rec
 real(dp), intent(in) :: mult,xmax
 real(dp),intent(inout) :: ene_out,ene_out1,ene_out2,ene_out3,ene_out4 !vz_i
!arrays
 real(dp), intent(in) :: an(0:nrec),bn2(0:nrec)
 
!Local variables-------------------------------
!scalars
 integer, parameter :: level = 7
 integer, save :: first = 1
 integer :: ii,kk,n_pt_integ_path2
 real(dp) :: epsilon,step,twotrotter,xmin,dr_step
 complex(dpc) :: D,Dnew,Dold,N,Nnew,Nold,dz_path,ene_acc,ene_acc1,ene_acc2
 complex(dpc) :: ene_acc3,ene_acc4
 complex(dpc) :: z_path,delta_calc
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 
! *************************************************************************

 
 call timab(611,1,tsec)
 
!structured debugging if debug_rec=T : print detailled result the first time we enter gran_potrec
 if(debug_rec .and. first==1)then
   write(message,'(a)')' ' 
   call wrtout(std_out,message,'PERS')
   write(message,'(a)')' gran_potrec : enter ' 
   call wrtout(std_out,message,'PERS')
   write(message,'(a,i8)')'n_pt_integ ' , n_pt_integ
   call wrtout(std_out,message,'COLL')
   first=0
 end if
 
 ene_out = zero
 ene_acc = czero
 ene_acc1 = czero
 ene_acc2 = czero
 ene_acc3 = czero
 ene_acc4 = czero
 

!path parameters 
!n_pt_integ = 2500
 xmin = -half
 step = (xmax-xmin)/real(n_pt_integ,dp)
 if(trotter==0)then
   twotrotter = one
   epsilon = .5d-1
 else
   twotrotter = two*real(trotter,dp)
   epsilon = half*sin( pi/twotrotter)
 end if

!xmin = -abs(xmin)**(1.d0/twotrotter)
 
!####################################################################
![xmax + i*epsilon,xmin + i*epsilon]
 dr_step = one/real(n_pt_integ,dp)
 dz_path = -cmplx((xmax-xmin)*dr_step,zero,dp)
 path1:  do ii = 0,n_pt_integ
!  z_path = cmplx(xmin + real(ii,dp)*(xmax-xmin)*dr_step,epsilon,dp)
   z_path = cmplx(xmin,epsilon,dp) - real(ii,dp)*dz_path
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
   delta_calc = func_rec(z_path,twotrotter)* N/D *dz_path
   if(ii==0.or.ii==n_pt_integ)then
     ene_acc = ene_acc + half*delta_calc
     if(debug_rec)  ene_acc1 = ene_acc1 + half*delta_calc
   else
     ene_acc = ene_acc + delta_calc
     if(debug_rec)  ene_acc1 = ene_acc1 + delta_calc
   end if
 end do path1
 
!####################################################################
![xmin + i*epsilon,xmin]
 if(epsilon/step>4.d0)then
   n_pt_integ_path2 = int(epsilon/step)+1
 else
   n_pt_integ_path2 = 5
 end if
 n_pt_integ_path2 = n_pt_integ
 dr_step = one/real(n_pt_integ_path2,dp)
 dz_path = -cmplx(zero,epsilon*dr_step,dp)
 path2:  do ii = 0,n_pt_integ_path2
!  z_path = cmplx(xmin,real(ii,dp)*epsilon*dr_step,dp)
   z_path = cmplx(xmin,zero,dp)-dz_path*real(ii,dp)
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
   delta_calc = func_rec(z_path,twotrotter)* N/D *dz_path             
   if(ii==0.or.ii==n_pt_integ_path2)then
     ene_acc = ene_acc + half*delta_calc           
     if(debug_rec) ene_acc3 = ene_acc3 + half*delta_calc
   else
     ene_acc = ene_acc + delta_calc
     if(debug_rec) ene_acc3 = ene_acc3 + delta_calc
   end if
 end do path2

 

!####################################################################
![xmin,0]
 if(xmin/=czero)then
   dr_step = one/real(n_pt_integ,dp)
   dz_path = cmplx(xmin*dr_step,zero,dp)
   path3:  do ii = 1,n_pt_integ !the integrand is 0 at 0
!    z_path = cmplx(real(ii,dp)*xmin*dr_step,zero,dp)
     z_path = real(ii,dp)*dz_path

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

!    <r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
     delta_calc = func_rec(z_path,twotrotter) * N/D *dz_path             
     if(ii==n_pt_integ)then
       ene_acc = ene_acc +half*delta_calc                                 
       if(debug_rec) ene_acc4 = ene_acc4 + half*delta_calc
     else
       ene_acc = ene_acc + delta_calc
       if(debug_rec) ene_acc4 = ene_acc4 +delta_calc
     end if
   end do path3
 end if
 
!####################################################################
![xmax,xmax+i*epsilon]
 dr_step = one/real(n_pt_integ_path2,dp)
 dz_path = cmplx(zero,epsilon*dr_step,dp)
 path4:  do ii = 0,n_pt_integ_path2
!  z_path = cmplx(xmax,real(ii,dp)*epsilon*dr_step,dp)
   z_path = cmplx(xmax,0,dp)+real(ii,dp)*dz_path

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
   delta_calc = func_rec(z_path,twotrotter) * N/D *dz_path               
   if(ii==0.or.ii==n_pt_integ_path2)then
     ene_acc = ene_acc + half*delta_calc                                
     if(debug_rec) ene_acc2 = ene_acc2 + half*delta_calc
   else
     ene_acc = ene_acc + delta_calc
     if(debug_rec) ene_acc2 = ene_acc2 + delta_calc
   end if
 end do path4
 
 ene_out = mult*real(ene_acc*cmplx(zero,-piinv,dp),dp)
 if(debug_rec) then
   ene_out1 = mult*real(ene_acc1*cmplx(zero,-piinv,dp),dp)
   ene_out2 = mult*real(ene_acc2*cmplx(zero,-piinv,dp),dp)
   ene_out3 = mult*real(ene_acc3*cmplx(zero,-piinv,dp),dp)
   ene_out4 = mult*real(ene_acc4*cmplx(zero,-piinv,dp),dp)
 end if

 call timab(611,2,tsec)

 contains 
 
!func_rec(z_path,twotrotter) = log(cone+z_path**twotrotter)

   function func_rec(z,x)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'func_rec'
!End of the abilint section

   implicit none

   complex(dpc) :: func_rec
   complex(dpc),intent(in) :: z
   real(dp),intent(in) :: x

   func_rec = log(cone+z**x)

 end function func_rec

end subroutine gran_potrec
!!***

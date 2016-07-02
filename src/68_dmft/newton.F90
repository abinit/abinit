!{\src2tex{textfont=tt}}
!!****f* ABINIT/newton
!! NAME
!! newton
!!
!! FUNCTION
!!  Compute root of a function with newton methods (newton/halley)
!! 
!! COPYRIGHT
!! Copyright (C) 2006-2016 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  x_input      : input of x             
!!  x_precision  : required precision on x
!!  max_iter     : maximum number of iterations
!!  opt_noninter
!!
!! OUTPUT
!!  f_precision  : output precision on function F
!!  ierr_hh      : different from zero if an error occurs
!!
!! PARENTS
!!      fermi_green
!!
!! CHILDREN
!!      compute_green,integrate_green
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"


subroutine newton(cryst_struc,green,paw_dmft,pawang,self,&
& x_input,x_precision,max_iter,f_precision,ierr_hh,opt_noninter,opt_algo)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_pawang,    only : pawang_type
 use m_crystal,   only : crystal_t
 use m_green,     only : green_type
 use m_paw_dmft,  only: paw_dmft_type
 use m_self,      only : self_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'newton'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 !type(MPI_type), intent(in) :: mpi_enreg
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(pawang_type),intent(in) :: pawang
 type(self_type), intent(inout) :: self
 integer,intent(in) :: opt_noninter,max_iter
 integer,intent(out) :: ierr_hh
 real(dp),intent(inout) :: x_input,x_precision
 real(dp),intent(inout) :: f_precision
 real(dp),intent(in), optional :: opt_algo
!Local variables-------------------------------
 integer iter
 real(dp) Fx,Fxprime,Fxdouble,xold,x_before,Fxoptimum
 real(dp) :: nb_elec_x
 integer option,optalgo
 logical l_minus,l_plus
 real(dp) :: x_minus,x_plus,x_optimum
 character(len=500) :: message
! *********************************************************************
 x_minus=-10_dp
 x_plus=-11_dp
 xold=-12_dp

 if(present(opt_algo)) then
   optalgo=opt_algo
 else
   optalgo=1 ! newton
 end if

 x_input=paw_dmft%fermie
 ierr_hh=0
 option =2  ! Halley method
 option =1  ! Newton method

!write(std_out,*) "ldaprint",opt_noninter

!--- Start of iterations
 write(message,'(a,3a)') "     Fermi level   Charge    Difference"
 call wrtout(std_out,message,'COLL')
!do iter=1, 40
!x_input=float(iter)/100_dp
!call function_and_deriv(cryst_struc,f_precision,green,iter,mpi_enreg,paw_dmft,pawang,self,&
!&  x_input,x_before,x_precision,Fx,Fxprime,Fxdouble,opt_noninter,option)
!write(std_out,*) x_input,Fx
!enddo
!call leave_new('COLL')

 l_minus=.false.
 l_plus=.false.
 Fxoptimum=1_dp
!========================================
!start iteration to find fermi level
!========================================
 do iter=1, max_iter

!  ========================================
!  If zero is located between two values: apply newton method or dichotomy
!  ========================================
   if(l_minus.and.l_plus) then

!    ==============================================
!    Compute the function and derivatives for newton
!    ==============================================
     call function_and_deriv(cryst_struc,f_precision,green,iter,paw_dmft,pawang,self,&
&     x_input,x_before,x_precision,Fx,Fxprime,Fxdouble,opt_noninter,option)

!    Apply stop criterion on Fx
     if(abs(Fx) < f_precision) then
!      write(message,'(a,2f12.6)') "Fx,f_precision",Fx,f_precision
!      call wrtout(std_out,message,'COLL')
       x_precision=x_input-xold
       return
     end if 
     if(iter==max_iter) then
       write(message,'(a,2f12.6)') "   Fermi level could not be found"
       call wrtout(std_out,message,'COLL')
       x_input=x_optimum
       ierr_hh=-123
       return
     end if

!    Cannot divide by Fxprime if too small
     if(abs(Fxprime) .le. 1.e-15)then
       ierr_hh=-314
       write(message,'(a,f12.7)') "Fxprime=",Fxprime
       call wrtout(std_out,message,'COLL')
       return
     end if

     x_precision=x_input-xold

!    ==============================================
!    Newton/Halley's  formula for next iteration
!    ==============================================
     xold=x_input
     if(option==1) x_input=x_input-Fx/Fxprime
     if(option==2) x_input=x_input-2*Fx*Fxprime/(2*Fxprime**2-Fx*Fxdouble)

!    ==============================================
!    If newton does not work well, use dichotomy.
!    ==============================================
     if(x_input<x_minus.or.x_input>x_plus) then
       call compute_nb_elec(cryst_struc,green,paw_dmft,pawang,self,&
&       Fx,opt_noninter,nb_elec_x,xold)
       write(message,'(a,3f12.6)') " ---",x_input,nb_elec_x,Fx
       call wrtout(std_out,message,'COLL')
       if(Fx>0) then
         x_plus=xold
       else if(Fx<0) then
         x_minus=xold
       end if
       x_input=(x_plus+x_minus)/2.d0
       
     end if
!    write(std_out,'(a,2f12.6)') " Q(xold) and dQ/dx=",Fx,Fxprime
!    write(std_out,'(a,f12.6)') " =>  new Fermi level",x_input
!    ========================================
!    Locate zero between two values
!    ========================================
   else 
     call compute_nb_elec(cryst_struc,green,paw_dmft,pawang,self,&
&     Fx,opt_noninter,nb_elec_x,x_input)
     write(message,'(a,3f12.6)') "  --",x_input,nb_elec_x,Fx
!    Possible improvement for large systems, removed temporarely for
!    automatic tests: more study is necessary: might worsen the convergency
!    if(iter==1) then
!    f_precision=max(abs(Fx/50),f_precision)
!    write(message,'(a,4x,a,e12.6)') ch10," Precision required changed to:",f_precision
!    call wrtout(std_out,message,'COLL')
!    endif
     call wrtout(std_out,message,'COLL')
     if(Fx>0) then
       l_plus=.true.
       x_plus=x_input
       x_input=x_input-0.02
     else if(Fx<0) then
       l_minus=.true.
       x_minus=x_input
       x_input=x_input+0.02
     end if
     
   end if

   if(abs(Fx)<abs(Fxoptimum)) then
     Fxoptimum=Fx
     x_optimum=x_input
   end if



!  if(myid==master) then
!  write(std_out,'(a,i4,3f12.6)') "i,xnew,F,Fprime",i,x_input,Fx,Fxprime
!  endif


!  ! Apply stop criterion on x
!  if(abs(x_input-xold) .le. x_input*x_precision) then
!  !    write(std_out,'(a,4f12.6)') "x_input-xold, x_precision*x_input   "&
!  !&    ,x_input-xold,x_precision*x_input,x_precision
!  f_precision=Fx
!  return
!  endif

 end do
!--- End of iterations


 ierr_hh=1
 return

 CONTAINS  !========================================================================================
!-----------------------------------------------------------------------
!!***

!!****f* newton/function_and_deriv
!! NAME
!!  function_and_deriv 
!!
!! FUNCTION
!!  Compute value of a function and its numerical derivatives
!!
!! INPUTS
!!  x_input      : input of x             
!!  option       : if 1 compute only first derivative
!!                 if 2 compute the two first derivatives.
!!  opt_noninter
!!
!! OUTPUTS
!!  Fx           : Value of F(x)
!!  Fxprime      : Value of F'(x)
!!  Fxdouble     : Value of F''(x)
!!
!! PARENTS
!!      newton
!!
!! CHILDREN
!!      compute_green,integrate_green
!!
!! SOURCE
subroutine function_and_deriv(cryst_struc,f_precision,green,iter,paw_dmft,pawang,self&
& ,x_input,x_old,x_precision,Fx,Fxprime,Fxdouble,opt_noninter,option)

 use m_profiling_abi

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_crystal, only : crystal_t
 use m_green, only : green_type
 use m_paw_dmft, only: paw_dmft_type
 use m_self, only : self_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'function_and_deriv'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 !type(MPI_type), intent(in) :: mpi_enreg
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(pawang_type),intent(in) :: pawang
 type(self_type), intent(inout) :: self
 integer,intent(in) :: iter,opt_noninter,option
 real(dp),intent(inout)  :: f_precision,x_input,x_precision
 real(dp),intent(out) :: Fx,Fxprime,Fxdouble
 real(dp),intent(inout) :: x_old
!Local variables-------------------------------
 real(dp) :: deltax,nb_elec_x,Fxmoins,Fxplus,xmoins,xplus,x0
 character(len=500) :: message
! *********************************************************************

!  choose deltax: for numeric evaluation of derivative
   if(iter==1) then
!    deltax=0.02
   end if
!  deltax=max((x_input-x_old)/10.d0,min(0.00001_dp,x_precision/100_dp))
   deltax=min(0.00001_dp,x_precision/100_dp)  ! small but efficient
!  endif
!  write(std_out,*) "iter,x_input,deltax",iter,x_input,deltax
   x0=x_input
   xmoins=x0-deltax
   xplus=x0+deltax

   call compute_nb_elec(cryst_struc,green,paw_dmft,pawang,self,&
&   Fx,opt_noninter,nb_elec_x,x0)

   write(message,'(a,3f12.6)') "  - ",x0,nb_elec_x,Fx
   call wrtout(std_out,message,'COLL')
!  write(std_out,*) "Fx", Fx
   if(abs(Fx)<f_precision) return

   call compute_nb_elec(cryst_struc,green,paw_dmft,pawang,self,&
&   Fxplus,opt_noninter,nb_elec_x,xplus)

   write(message,'(a,3f12.6)') "  - ",xplus,nb_elec_x,Fxplus
   call wrtout(std_out,message,'COLL')

   if(option==2) then
     call compute_nb_elec(cryst_struc,green,paw_dmft,pawang,self,&
&     Fxmoins,opt_noninter,nb_elec_x,xmoins)

     write(message,'(a,3f12.6)') "  - ",xmoins,nb_elec_x,Fxmoins
     call wrtout(std_out,message,'COLL')
   end if

   if(option==1) then
     Fxprime=(Fxplus-Fx)/deltax
   else if (option==2) then
     Fxprime=(Fxplus-Fxmoins)/(2*deltax)
     Fxdouble=(Fxplus+Fxmoins-2*Fx)/(deltax**2)
   end if
!  write(std_out,*) "after computation of Fxprime",myid
   if(Fxprime<zero) then
     write(message,'(a,f12.6)') "  Warning: slope of charge versus fermi level is negative !", Fxprime
     call wrtout(std_out,message,'COLL')
   end if
   x_old=x_input

   return
 end subroutine function_and_deriv
!!***

!!****f* newton/compute_nb_elec
!! NAME
!! compute_nb_elec
!!
!! FUNCTION
!!  Compute nb of electrons as a function of Fermi level
!!
!! INPUTS
!! fermie       : input of energy 
!! opt_noninter
!!
!! OUTPUTS
!! Fx           : Value of F(x).
!! nb_elec_x    : Number of electrons for the value of x
!!
!! PARENTS
!!      newton
!!
!! CHILDREN
!!      compute_green,integrate_green
!!
!! SOURCE

subroutine compute_nb_elec(cryst_struc,green,paw_dmft,pawang,self,&
&  Fx,opt_noninter,nb_elec_x,fermie)

 use m_profiling_abi

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_crystal, only : crystal_t
 use m_green, only : green_type,compute_green,integrate_green
 use m_paw_dmft, only: paw_dmft_type
 use m_self, only : self_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'compute_nb_elec'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 !type(MPI_type), intent(in) :: mpi_enreg
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(pawang_type),intent(in) :: pawang
 type(self_type), intent(inout) :: self
 integer,intent(in) :: opt_noninter
 real(dp),intent(in)  :: fermie
 real(dp),intent(out) :: Fx,nb_elec_x
! *********************************************************************
   paw_dmft%fermie=fermie
   call compute_green(cryst_struc,green,paw_dmft,pawang,0,self,opt_self=1,&
&   opt_nonxsum=1,opt_nonxsum2=1)
   call integrate_green(cryst_struc,green,paw_dmft,pawang,prtopt=0,&
&   opt_ksloc=-1) !,opt_nonxsum=1)
!  opt_ksloc=-1, compute total charge
   nb_elec_x=green%charge_ks
   Fx=nb_elec_x-paw_dmft%nelectval

   if(opt_noninter==1) then
   end if
 end subroutine compute_nb_elec
!!***

end subroutine newton
!!***

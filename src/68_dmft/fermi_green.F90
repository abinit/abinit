!{\src2tex{textfont=tt}}
!!****f* ABINIT/fermi_green
!! NAME
!! fermi_green
!!
!! FUNCTION
!!  Compute Fermi level for DMFT or LDA.
!! 
!! COPYRIGHT
!! Copyright (C) 2006-2016 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! fermie: current value
!! f_precision: precision of f required
!! ITRLST: =1 if last iteration of DMFT
!! opt_noninter   : if one wants the LDA fermi level
!! max_iter : max number of iterations.
!!
!! OUTPUT
!! fermie: output value
!!
!! CHILDREN
!!  wrtout, spline_c,invfourier, nfourier
!!
!! PARENTS
!!      dmft_solve
!!
!! CHILDREN
!!      compute_green,integrate_green,newton,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

subroutine fermi_green(cryst_struc,green,paw_dmft,pawang,&
& self)

 use m_profiling_abi

 use defs_basis
 use defs_abitypes
 use m_pawang, only : pawang_type
 use m_crystal, only : crystal_t
 use m_green, only : green_type,integrate_green,compute_green
 use m_paw_dmft, only: paw_dmft_type
 use m_self, only : self_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fermi_green'
 use interfaces_14_hidewrite
 use interfaces_68_dmft, except_this_one => fermi_green
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 type(paw_dmft_type), intent(inout) :: paw_dmft
 !type(MPI_type), intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(self_type), intent(inout) :: self

!Local variables-------------------------------
 integer :: ierr_hh,opt_noninter,max_iter
 real(dp):: x_precision,f_precision,fermi_old
! real(dp) :: hx 
 character(len=500) :: message
!************************************************************************
!
 write(message,'(a,8x,a)') ch10,"  == Compute Fermi level"
 call wrtout(std_out,message,'COLL')

!=============
!headers
!=============
 write(message,'(2a)') ch10, "  |---Newton method to search Fermi level ------------|"
 call wrtout(std_out,message,'COLL')
 write(message,'(2a,f13.6)') ch10, "  |--- Initial value for Fermi level",paw_dmft%fermie
 call wrtout(std_out,message,'COLL')

!========================================
!Define precision and nb of iterations
!=========================================
 fermi_old=paw_dmft%fermie
 ierr_hh=0
 f_precision=paw_dmft%dmft_chpr
 !f_precision=0.01
 x_precision=tol5
!if(option==1) then
!f_precision=(erreursurlacharge)/100d0
!else
!f_precision=tol11
!endif
 max_iter=1 ! for tests only
 !write(6,*) "for tests max_iter=1"
 max_iter=50
 opt_noninter=4

!=====================
!Call newton method
!=====================
 write(message,'(a,4x,a,e13.6)') ch10," Precision required :",f_precision
 call wrtout(std_out,message,'COLL')
 if (f_precision<10_dp)  then
   call newton(cryst_struc,green,paw_dmft,pawang,self,&
&   paw_dmft%fermie,x_precision,max_iter,f_precision,ierr_hh,opt_noninter)
 end if

!===========================
!Deals with errors signals
!===========================
 if(ierr_hh==-314) then
   write(message,'(a)') "Warning, check Fermi level"
   call wrtout(std_out,message,'COLL')
!  call leave_new('COLL')
   write(message,'(2a,f13.6)') ch10, "  |---  Final  value for Fermi level (check)",paw_dmft%fermie
   call wrtout(std_out,message,'COLL')
 else if (ierr_hh==-123) then
   write(message,'(a,f13.6)') " Fermi level is put to",fermi_old
   paw_dmft%fermie=fermi_old
   call wrtout(std_out,message,'COLL')

!  =====================================
!  If fermi level search was successful
!  =====================================
 else
   write(message,'(a,4x,a,e13.6)') ch10," Precision achieved on Fermi Level :",x_precision
   call wrtout(std_out,message,'COLL')
   write(message,'(4x,a,e13.6)') " Precision achieved on number of electrons :",f_precision
   call wrtout(std_out,message,'COLL')
   write(message,'(2a,f13.6)') ch10, "  |---  Final  value for Fermi level",paw_dmft%fermie
   call wrtout(std_out,message,'COLL')
 end if

!========================================================
!Check convergence of fermi level during DMFT iterations
!========================================================
 if(paw_dmft%idmftloop>=2) then
   if(abs(paw_dmft%fermie-fermi_old).le.paw_dmft%dmft_fepr) then
!    write(message,'(a,8x,a,e9.2,a,8x,a,e12.5)') ch10,"|fermie(n)-fermie(n-1)|=<",paw_dmft%dmft_fepr,ch10,&
     write(message,'(a,8x,a,e9.2,a,e9.2,a,8x,a,e12.5)') ch10,"|fermie(n)-fermie(n-1)|=",&
&     abs(paw_dmft%fermie-fermi_old),"<",paw_dmft%dmft_fepr,ch10,&
&     "=> DMFT Loop: Fermi level is converged to:",paw_dmft%fermie
     call wrtout(std_out,message,'COLL')
     green%ifermie_cv=1
   else
     write(message,'(a,8x,a,2f12.5)') ch10,"DMFT Loop: Fermi level is not converged:",&
&     paw_dmft%fermie
     call wrtout(std_out,message,'COLL')
     green%ifermie_cv=0
   end if
 end if
 write(message,'(2a)') ch10, "  |---------------------------------------------------|"
 call wrtout(std_out,message,'COLL')
!

!==========================================================
!Recompute full green function including non diag elements
!==========================================================
 call compute_green(cryst_struc,green,paw_dmft,pawang,0,self,opt_self=1,opt_nonxsum=1)
 call integrate_green(cryst_struc,green,paw_dmft,pawang,prtopt=0,opt_ksloc=3) !,opt_nonxsum=1)

 return 
end subroutine fermi_green
!!***

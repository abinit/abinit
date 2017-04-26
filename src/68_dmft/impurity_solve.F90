!{\src2tex{textfont=tt}}
!!****f* ABINIT/impurity_solve
!! NAME
!! impurity_solve
!!
!! FUNCTION
!! Solve the Impurity problem
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  lda_occup
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawtab <type(pawtab)>
!!
!! OUTPUT
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!
!! NOTES
!!
!! PARENTS
!!      dmft_solve
!!
!! CHILDREN
!!      copy_green,destroy_green_tau,flush_unit,fourier_green,hubbard_one
!!      init_green_tau,integrate_green,ldau_self,print_green,print_matlu
!!      printocc_green,qmc_prep_ctqmc,timab,trace_oper,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

subroutine impurity_solve(cryst_struc,green,hu,paw_dmft,&
& pawang,pawtab,self_old,self_new,weiss,pawprtvol)


 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi

 use m_crystal, only : crystal_t
 use m_green, only : green_type, fourier_green&
& ,init_green_tau,destroy_green_tau,print_green,printocc_green,integrate_green,copy_green
 use m_paw_dmft, only : paw_dmft_type
 use m_oper,     only : trace_oper
 use m_hu, only : hu_type
 use m_matlu, only : matlu_type,print_matlu,init_matlu,destroy_matlu
 use m_self, only : self_type
 use m_energy, only : energy_type
 use m_pawang, only : pawang_type
 use m_pawtab, only : pawtab_type
 use m_io_tools, only : flush_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'impurity_solve'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_68_dmft, except_this_one => impurity_solve
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
! type(pawang_type), intent(in) :: pawang
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type), intent(inout) :: weiss
 type(green_type), intent(inout) :: green !vz_i
 type(hu_type),intent(inout) :: hu(cryst_struc%ntypat)
 !type(MPI_type), intent(in) :: mpi_enreg
 type(pawang_type), intent(in) :: pawang
 type(pawtab_type),intent(in)  :: pawtab(cryst_struc%ntypat)
 type(paw_dmft_type), intent(inout)  :: paw_dmft
 type(self_type), intent(inout) :: self_new
 type(self_type), intent(inout) :: self_old
 integer, intent(in) :: pawprtvol

!Local variables ------------------------------
 real(dp) :: tsec(2)
 character(len=500) :: message
 complex(dpc) :: xx
 integer :: ifreq   
! integer iatom,il,i_nd,isppol,lpawu,im,Nd,nrat,nsweeptot
! real(dp) :: acc,kx
! real(dp), allocatable :: correl(:,:),g0(:,:),gtmp(:,:)
!scalars
!************************************************************************
!character(len=500) :: message

 call timab(622,1,tsec)
!=======================================================================
!== Prepare data for Hirsch Fye QMC
!== NB: for CTQMC, Fourier Transformation are done inside the CTQMC code
!=======================================================================
 if(abs(paw_dmft%dmft_solv)==4) then
!  == Initialize weiss and green functions for fourier transformation
!  -------------------------------------------------------------------
   write(message,'(2a,i3,13x,a)') ch10,'   ===  Initialize Weiss field G_0(tau)'
   call wrtout(std_out,message,'COLL')
   call init_green_tau(weiss,paw_dmft)
   call init_green_tau(green,paw_dmft)
!  in init_solver

!  == Print weiss function G_0(tau=0-) before computation (really useless check)
!  ------------------------------------------------------------------------------
   if(abs(pawprtvol)>3) then
     write(message,'(2a,i3,13x,a)') ch10,'   ===  Check G_0(tau=0-) first'
     call wrtout(std_out,message,'COLL')
     call printocc_green(weiss,6,paw_dmft,3)
   end if

!  == Fourier transform of weiss Field
!  ------------------------------------
!  for fourier of KS green functions
!  call fourier_green(cryst_struc,weiss,mpi_enreg,paw_dmft,pawang,pawtab,1)
   write(message,'(2a,i3,13x,a)') ch10,'   ===  Inverse Fourier Transform w->t of Weiss Field'
   call wrtout(std_out,message,'COLL')
   call fourier_green(cryst_struc,weiss,paw_dmft,pawang,opt_ksloc=2,opt_tw=-1)

!  == Print weiss function G2_0(tau=0-) 
!  --------------------------------------
   call printocc_green(weiss,6,paw_dmft,3,opt_weissgreen=1)

!  for fourier of KS green functions
!  call fourier_green(cryst_struc,weiss,mpi_enreg,paw_dmft,pawang,pawtab,1)
!  == Print G_0(tau) in files
!  ---------------------------
   if(paw_dmft%dmft_prgn==1) then
     call print_green('weiss',weiss,1,paw_dmft,pawprtvol=1,opt_wt=2)
   end if

 else if(abs(paw_dmft%dmft_solv)>=5) then
!  == Initialize  green functions for imaginary times
!  -------------------------------------------------------------------
   write(message,'(2a,i3,13x,a)') ch10,'   ===  Initialize Green function G(tau)'
   call wrtout(std_out,message,'COLL')
   call init_green_tau(green,paw_dmft)
   
 end if
!=======================================================================
!== End preparation of QMC
!=======================================================================

!=======================================================================
!== Solve impurity model   =============================================
!=======================================================================
 write(message,'(2a,i3,13x,a)') ch10,'  ===  Solve impurity model'
 call wrtout(std_out,message,'COLL')
 if(abs(paw_dmft%dmft_solv)==1) then

!  == LDA+U for test -> self
!  -------------------
   call ldau_self(cryst_struc,green,paw_dmft,&
&   pawtab,self_new,opt_ldau=1,prtopt=pawprtvol)
 else if(abs(paw_dmft%dmft_solv)==2) then

!  == Hubbard One -> green
!  -------------------
   call hubbard_one(cryst_struc,green,hu,paw_dmft,&
&   pawang,pawprtvol,self_old%hdc,weiss)

 else if(abs(paw_dmft%dmft_solv)==4) then

!  == QMC
!  -------------------
   call copy_green(weiss,green,opt_tw=1)
!  call qmc_prep
   message = '  ===  QMC not yet distributed '
   MSG_ERROR(message)
!   call qmc_prep(cryst_struc,green,hu,mpi_enreg,paw_dmft,pawang&
!&   ,pawprtvol,self_old%qmc_xmu,weiss)

 else if(abs(paw_dmft%dmft_solv)>=5) then

!  == Nothing
!  -------------------
!   call copy_green(weiss,green,opt_tw=1)
!   call copy_green(weiss,green,opt_tw=2)

   call qmc_prep_ctqmc(cryst_struc,green,self_old,hu,paw_dmft,pawang,pawprtvol,weiss)


 else if(abs(paw_dmft%dmft_solv)==0) then

!  == Nothing
!  -------------------
!  weiss%occup%has_operks=0 -> only local part is duplicated
   call copy_green(weiss,green,opt_tw=2)
 end if
!call print_green("invWeiss",cryst_struc,weiss,3,paw_dmft,pawtab,2)

!=======================================================================
!== Treat data from HF QMC
!=======================================================================
 if(abs(paw_dmft%dmft_solv)>=4) then
!  propagate qmc_shift (useful for compute_energy)
   if(abs(paw_dmft%dmft_solv)==4) then
     self_new%qmc_shift(:)=self_old%qmc_shift(:)
     self_new%qmc_xmu(:)=self_old%qmc_xmu(:)
   end if

!  == Print local occupations from G(tau)
!  ---------------------------------------

!  == Fourier back transform of green function G(tau)->G(iw_n) and
!  == compute occupations from g(tau)
!  -------------------------------------------------------------------
   if(abs(paw_dmft%dmft_solv)==4) then
     write(message,'(2a,i3,13x,a)') ch10,'   ===  Direct Fourier Transform t->w of Green Function'
     call wrtout(std_out,message,'COLL')
     call fourier_green(cryst_struc,green,paw_dmft,&
&     pawang,opt_ksloc=2,opt_tw=1)
     do ifreq=1,green%nw
       xx= green%oper(ifreq)%matlu(1)%mat(1,1,1,1,1) 
       write(112,*) paw_dmft%omega_lo(ifreq),real(one/xx),aimag(one/xx)
       write(113,*) paw_dmft%omega_lo(ifreq),real(xx),aimag(xx)
     end do
     call flush_unit(112)
     call flush_unit(113)
!   if(paw_dmft%dmft_solv==5) stop
     if(pawprtvol>=3) then
       write(message,'(a,2x,a,f13.5)') ch10,&    ! debug
&      " == Print green function for small freq after fourier " ! debug
       call wrtout(std_out,message,'COLL')    ! debug
       call print_matlu(green%oper(1)%matlu,paw_dmft%natom,1)    ! debug
     end if

     write(message,'(2a,i3,13x,a)') ch10,'   INVERSE FOURIER OF G0 SUPPRESSED'
     call wrtout(std_out,message,'COLL')
   end if
   if(abs(paw_dmft%dmft_solv)==888) then
!  == Back fourier transform of G_0(tau) for compensation (try or comment or improve FT).
!  -------------------------------------------------------------------
     write(message,'(2a,i3,13x,a)') ch10,'   ===  Direct Fourier transform t->w of Weiss'
     call wrtout(std_out,message,'COLL')
     call fourier_green(cryst_struc,weiss,paw_dmft,&
&     pawang,opt_ksloc=2,opt_tw=1)

     if(pawprtvol>=3) then
       write(message,'(a,2x,a,f13.5)') ch10,&    ! debug
&      " == Print weiss function for small freq after fourier " ! debug
       call wrtout(std_out,message,'COLL')    ! debug
       call print_matlu(weiss%oper(1)%matlu,paw_dmft%natom,1)    ! debug
     end if
     call destroy_green_tau(weiss)
   end if

!  == Destroy tau part of green
!  -------------------------------------------------------------------
   call trace_oper(green%occup_tau,green%charge_ks,green%charge_matlu_solver,2)
   green%has_charge_matlu_solver=2
   call destroy_green_tau(green)
 end if
!=======================================================================
!== End Treat data for QMC
!=======================================================================

!=======================================================================
!== Integrate green function and printout occupations
!=======================================================================
!For dmft_solv=-1,0,or 1 , the green function was not yet computed: it
!cannot be integrated
!=======================================================================
 if(paw_dmft%dmft_solv>=2.and.green%w_type=="imag") then
!  ==  Integrate G(iw_n)
!  ---------------------
   write(message,'(2a,i3,13x,a)') ch10,'   ===  Integrate local part of green function'
   call wrtout(std_out,message,'COLL')
   call integrate_green(cryst_struc,green,paw_dmft,&
&   pawang,prtopt=2,opt_ksloc=2,opt_after_solver=1)

!  == Print local occupations from integration of G(iw_n)
!  --------------------------------------------------------
   call printocc_green(green,5,paw_dmft,3)

!  == Print G_loc(w)
!  --------------------------------------------------------
   if(paw_dmft%dmft_prgn==1) then
     call print_green('DMFT_IMPURITY',green,1,paw_dmft,pawprtvol=1,opt_wt=1)
   end if
 end if
!stop 

 if(abs(pawprtvol)>0) then
 end if

 call timab(622,2,tsec)
end subroutine impurity_solve
!!***

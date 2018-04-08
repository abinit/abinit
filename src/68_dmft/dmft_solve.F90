!{\src2tex{textfont=tt}}
!!****f* ABINIT/dmft_solve
!! NAME
!! dmft_solve
!!
!! FUNCTION
!! Solve the DMFT loop from PAW data.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  istep           =  step of iteration for LDA.
!!  lda_occup <type(oper_type)> = occupations in the correlated orbitals in LDA
!!  paw_dmft <type(paw_dmft_type)> =  data for self-consistent LDA+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pawprtvol  = option for printing
!!
!! OUTPUT
!!  paw_dmft <type(paw_dmft_type)> =  data for self-consistent LDA+DMFT calculations.
!!
!! NOTES
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      check_fourier_green,compute_energy,compute_green,compute_ldau_energy
!!      data4entropydmft_setdc,data4entropydmft_sethu,dc_self,destroy_energy
!!      destroy_green,destroy_hu,destroy_self,diff_oper,dyson,fermi_green
!!      icip_green,impurity_solve,init_energy,init_green,init_hu
!!      initialize_self,integrate_green,local_ks_green,make_qmcshift_self
!!      new_self,print_self,printocc_green,psichi_renormalization,rw_self
!!      spectral_function,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

subroutine dmft_solve(cryst_struc,istep,lda_occup,paw_dmft,pawang,pawtab,pawprtvol)


 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_profiling_abi
 use m_data4entropyDMFT

 use m_time,           only : timab
 use m_pawang, only : pawang_type
 use m_pawtab, only : pawtab_type
 use m_paw_dmft, only: paw_dmft_type
 use m_crystal, only : crystal_t
 use m_green, only : green_type, destroy_green, icip_green,init_green,&
&                    print_green,printocc_green,&
&                    integrate_green,copy_green,compute_green,check_fourier_green
 use m_oper, only : oper_type,diff_oper,upfold_oper,loc_oper
 use m_self, only : self_type,initialize_self,destroy_self,print_self,dc_self,rw_self,new_self,make_qmcshift_self
 use m_hu, only : hu_type,init_hu,destroy_hu
 use m_energy, only : energy_type,init_energy,destroy_energy,compute_energy,print_energy,compute_ldau_energy
 use m_matlu, only : print_matlu,sym_matlu

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dmft_solve'
 use interfaces_14_hidewrite
 use interfaces_68_dmft, except_this_one => dmft_solve
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: istep
 integer, intent(in) :: pawprtvol
 !type(MPI_type), intent(inout) :: mpi_enreg
 type(pawang_type), intent(in) :: pawang
 type(crystal_t),intent(in) :: cryst_struc
 type(pawtab_type),intent(in)  :: pawtab(cryst_struc%ntypat)
 type(oper_type),intent(in)  :: lda_occup
 type(paw_dmft_type), intent(inout)  :: paw_dmft

!Local variables ------------------------------
!array
 real(dp) :: tsec(2)
!scalars
 integer :: check,idmftloop,istep_iter,spaceComm,my_rank,opt_renorm
 integer :: itypat
 logical :: etot_var
 character(len=200) :: char_enddmft
! type
 type(green_type) :: green
 type(green_type) :: greenlda
 type(hu_type),allocatable :: hu(:)
 type(green_type) :: weiss
 type(self_type) :: self
 type(self_type) :: self_new
 type(energy_type) :: energies_dmft
 type(energy_type) :: energies_tmp
 character(len=500) :: message
 character(len=5) :: thdyn
 character(len=4) :: part2,part3
!************************************************************************

 DBG_ENTER('COLL')
 my_rank = xmpi_comm_rank(paw_dmft%spacecomm)

 check=paw_dmft%dmftcheck ! checks enabled
 paw_dmft%dmft_fepr=tol5
 paw_dmft%dmft_chpr=tol6
!paw_dmft%dmft_chpr=20_dp ! total number of electron.
 paw_dmft%dmft_prgn=1
 paw_dmft%dmft_prgn=0
 etot_var=.true.
 thdyn="fcalc"
 thdyn="ecalc"
 if(thdyn=="ecalc") then ! valid
   part2="both"
   part3="none"
 else if(thdyn=="fcalc") then ! not tested
   part2="corr"
   part3="band"
 end if

 if(check==1) then
   write(message,'(2a)') ch10,' DMFT Checks are enabled '
 else
   write(message,'(2a)') ch10,' DMFT Checks will not be performed'
 end if
 call wrtout(std_out,message,'COLL')


 if(istep==0) then
   message = ' istep should not be equal to zero'
   MSG_BUG(message)
 end if
 spaceComm=paw_dmft%spacecomm
 !if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_kpt
 !call xmpi_barrier(spaceComm)

 call initialize_self(self,paw_dmft)
 call init_energy(cryst_struc,energies_dmft)

!===========================================================================
!==  First construct LDA green function (Init, Compute, Integrate, Print)
!===========================================================================
 write(message,'(6a)') ch10,' =====================================', &
& ch10,' =====  LDA Green Function Calculation',&
& ch10,' ====================================='
 call wrtout(std_out,message,'COLL')
 call icip_green("LDA",cryst_struc,greenlda,paw_dmft,pawang,3,self)
 !call print_green('LDA_NOT_renormalized',greenlda,1,paw_dmft,pawprtvol=1,opt_wt=1)

!== Compare greenlda%occup and lda_occup: check that LDA green function is fine
!----------------------------------------------------------------------
 write(message,'(2a)') ch10,&
& '  == Check lda occ. mat. from green with respect to the direct calc =='
 call wrtout(std_out,message,'COLL')
 call diff_oper("Occup from LDA green function",&
& "LDA occupations",greenlda%occup,lda_occup,1,paw_dmft%dmft_tolfreq)
! write(message,'(2a)') ch10,&
!& '  ***** => Warning : diff_oper is suppressed for test'
! call wrtout(std_out,message,'COLL')
 write(message,'(2a)') ch10,&
& '  ***** => Calculation of Green function is thus correct without self ****'
 call wrtout(std_out,message,'COLL')
 call destroy_green(greenlda)
 
!== Orthonormalize psichi
!----------------------------------------------------------------------
 call timab(621,1,tsec)
 opt_renorm=3
 if(paw_dmft%nspinor==2.and.paw_dmft%dmft_solv==5) opt_renorm=2 ! necessary to use hybri_limit in qmc_prep_ctqmc
                                                                ! ought to be  generalized  in the  future
 if(paw_dmft%dmft_solv/=-1) then
   call psichi_renormalization(cryst_struc,paw_dmft,pawang,opt=opt_renorm)

!  ===========================================================================
!  ==  re-construct LDA green function with new psichis
!  ===========================================================================
   write(message,'(6a)')&
&   ch10,' ==============================================================='&
&   ,ch10,' =====  LDA Green Function Calculation with renormalized psichi'&
&   ,ch10,' =============================================================='
 end if
 call timab(621,2,tsec)

 call wrtout(std_out,message,'COLL')
 call icip_green("LDA_renormalized",cryst_struc,greenlda,&
& paw_dmft,pawang,pawprtvol,self)
 !call print_green('LDA_renormalized',greenlda,1,paw_dmft,pawprtvol=1,opt_wt=1)

 !Need to store idmftloop and set it to zero to avoid useless print_energy in ab_out
 idmftloop=paw_dmft%idmftloop
 paw_dmft%idmftloop=0
 call compute_energy(cryst_struc,energies_dmft,greenlda,paw_dmft,pawprtvol,pawtab,self,occ_type=" lda",part='both')
 paw_dmft%idmftloop=idmftloop

 if(paw_dmft%dmft_prgn==1) then
   if(paw_dmft%lpsichiortho==1) then
     call local_ks_green(greenlda,paw_dmft,prtopt=1)
   end if
 end if
 call destroy_self(self)
!call printocc_green(greenlda,9,paw_dmft,3,chtype="LDA GREEN PSICHI")

 write(message,'(6a)')&
& ch10,' ==============================================================='&
& ,ch10,' ===== Define Interaction and self-energy' &
& ,ch10,' =============================================================='
 call wrtout(std_out,message,'COLL')
!== define Interaction from input upawu and jpawu
!----------------------------------------------------------------------
 ABI_DATATYPE_ALLOCATE(hu,(cryst_struc%ntypat))
 call init_hu(cryst_struc,pawtab,hu,paw_dmft%dmftqmc_t2g)
 call initialize_self(self,paw_dmft)

 ! Set Hu in density representation for calculation of entropy if needed...
 do itypat=1,cryst_struc%ntypat
   if ( hu(itypat)%lpawu /= -1 ) then
     call data4entropyDMFT_setHu(paw_dmft%forentropyDMFT,itypat,hu(itypat)%udens(:,:))
   end if
 end do

!== define self from scratch or file and double counting
!----------------------------------------------------------------------
!-  Self allocated
 call dc_self(greenlda%charge_matlu,cryst_struc,hu,self,paw_dmft%dmft_dc,pawprtvol)

!-   Read self or do self=hdc
 if(paw_dmft%dmft_solv==4) then
!  write(std_out,*) "shift before rw_self",self%qmc_shift(1)
   call make_qmcshift_self(cryst_struc,hu,self)
 end if
 call timab(627,1,tsec)
 call rw_self(self,paw_dmft,prtopt=2,opt_rw=1,istep_iter=1000*istep)
 call timab(627,2,tsec)

!== If QMC is used,  and self energy is read for file, then
!== one does NOT shifts the self-energy because it was already shifted when writed,
!==  and thus then weiss will be shifted
!----------------------------------------------------------------------
!if(paw_dmft%dmft_solv==4.and.paw_dmft%dmft_rslf==1) &
!&           call make_qmcshift_self(cryst_struc,hu,self)
!if(paw_dmft%dmft_solv==4.and.paw_dmft%dmft_rslf/=1) &
!&           call make_qmcshift_self(cryst_struc,hu,self,apply=.true.)

 call destroy_green(greenlda)  ! destroy LDA green function
 call print_self(self,"print_dc",paw_dmft,prtopt=2)



!===========================================================================
!==  construct green function with the self-energy.
!===========================================================================
 write(message,'(6a)') &
& ch10,' =================================================================', &
& ch10,' =====  Green Function Calculation with input self-energy ========', &
& ch10,' ================================================================='
 call wrtout(std_out,message,'COLL')
 call icip_green("Green_inputself",cryst_struc,green,&
& paw_dmft,pawang,pawprtvol,self,opt_self=1)
   !call print_green('beforefermi_green',green,1,paw_dmft,pawprtvol=1,opt_wt=1)
!   stop
!   call leave_new('COLL')

!== Find fermi level
!---------------------------------------------------------------------
!write(message,'(2a,i3,13x,a)') ch10,'   ===  Compute green function from self-energy'
 call fermi_green(cryst_struc,green,paw_dmft,pawang,self)
!== define weiss field only for the local quantities (opt_oper=2)
!----------------------------------------------------------------------
! write(std_out,*) "nkpt  befreo init_greenweiss",ifreq,paw_dmft%nkpt
 call init_green(weiss,paw_dmft,opt_oper_ksloc=2)
! do ifreq=1,weiss%nw
!   write(std_out,*) "nkpt from weiss1",ifreq,weiss%oper(ifreq)%nkpt
! enddo

!== Check fourier transforms
!----------------------------------------------------------------------
 if(check==1) then
   call check_fourier_green(cryst_struc,green,paw_dmft,pawang)
 end if

!== If QMC is used,  and self energy is not read for file, then
!== one shifts the self-energy, and thus then weiss will be shifted
!== after dyson, in a coherent way concerning qmc_shift and qmc_xmu.
!----------------------------------------------------------------------
!if(paw_dmft%dmft_solv==4.and.paw_dmft%dmft_rslf/=1) &
!&           call make_qmcshift_self(cryst_struc,hu,self,apply=.true.)
!if(paw_dmft%dmft_solv==4) write(std_out,*) "shift after make_qmcshift_self",self%qmc_shift(1)

 write(message,'(6a)') &
& ch10,' ======================================================'&
& ,ch10,' =====  DMFT Loop starts here                  ========'&
& ,ch10,' ======================================================'
 call wrtout(std_out,message,'COLL')

!=======================================================================
!===  dmft loop  =======================================================
 do idmftloop=1, paw_dmft%dmft_iter 
   !paw_dmft%idmftloop=idmftloop
   paw_dmft%idmftloop=paw_dmft%idmftloop+1
!  =======================================================================
   istep_iter=1000*istep+idmftloop

   write(message,'(2a,i3,13x,a)') ch10,&
&   ' =====  DMFT Loop : ITER number',paw_dmft%idmftloop,'========'
   call wrtout(std_out,message,'COLL')

!  == Dyson Equation G,self -> weiss(w)
!  ---------------------------------------------------------------------
   call dyson(green,paw_dmft,self,weiss,opt_weissself=1)  
!   call print_green('afterDyson',green,1,paw_dmft,pawprtvol=1,opt_wt=1)
!   stop
!   call leave_new('COLL')

!  == Printout local "occupations" from weiss field  (useless)
   if(abs(pawprtvol)>3) then
     call integrate_green(cryst_struc,weiss,paw_dmft,&
&     pawang,prtopt=2,opt_ksloc=2)
     call printocc_green(weiss,5,paw_dmft,3,opt_weissgreen=1)
   end if
   
!  ===  Prepare data, solve Impurity problem: weiss(w) -> G(w)
!  ---------------------------------------------------------------------
   call initialize_self(self_new,paw_dmft)

   call impurity_solve(cryst_struc,green,hu,&
&   paw_dmft,pawang,pawtab,self,self_new,weiss,pawprtvol) ! weiss-> green, or self if dmft_solv=1
!  if(paw_dmft%dmft_solv==4)  write(std_out,*) "shift after impurity",self%qmc_shift(1)

!  ==  Compute double counting from charge from green_solver
!  ---------------------------------------------------------------------
   if (green%has_charge_matlu_solver/=2) green%charge_matlu_solver=green%charge_matlu
   call dc_self(green%charge_matlu_solver,cryst_struc,hu,self_new,paw_dmft%dmft_dc,pawprtvol)

!  ==  Solve dyson equation. G_imp(w), weiss_imp(w) -> Self_imp(w)
!  ---------------------------------------------------------------------
!  if dmft_solv==1, self is computed previously
   if(abs(paw_dmft%dmft_solv)/=1) then
     call dyson(green,paw_dmft,self_new,weiss,opt_weissself=2) 
   end if
!  do ifreq=1,green%nw
!  call sym_matlu(cryst_struc,self%oper(ifreq)%matlu,pawang)
!  enddo

!  ==  Possibility if imposing self (opt_rw==3)
!  ---------------------------------------------------------------------
   call timab(627,1,tsec)
   call rw_self(self_new,paw_dmft,prtopt=2,opt_rw=3,istep_iter=istep_iter)
   call timab(627,2,tsec)

!  print dc just computed before and self computed in before in dyson or
!  impurity_solve
   if(abs(pawprtvol)>=3) then
     write(message,'(2a)') ch10,"  == New self"
     call wrtout(std_out,message,'COLL')
     call print_self(self_new,"print_dc",paw_dmft,2)
     write(message,'(2a)') ch10,"  == Old self"
     call wrtout(std_out,message,'COLL')
     call print_self(self,"print_dc",paw_dmft,2)
   end if

!  if(paw_dmft%dmft_solv==4) write(std_out,*) "shift before computeenergy ",self%qmc_shift(1)
!  ==  Compute Energy with NEW self-energy and edc from green_solver,
!  new local green function and old occupations for eband
!  fermi level not optimized for this self_energy.
!  ---------------------------------------------------------------------
!  green= local green function and local charge comes directly from solver
!  green= ks green function and occupations comes from old_self 
   call compute_energy(cryst_struc,energies_dmft,green,paw_dmft,pawprtvol,&
&   pawtab,self_new,occ_type="nlda",part=part2)

!  ==  Mix new and old self_energies and double countings
!  ---------------------------------------------------------------------
   call new_self(self,self_new,paw_dmft,1) ! self,self_new => self
   write(message,'(2a)') ch10,"  == After mixing,"
     !print *, " my_rank newself", my_rank,self%oper(1)%matlu(1)%mat(1,1,1,1,1)
   call wrtout(std_out,message,'COLL')
   call print_self(self,"print_dc",paw_dmft,2) ! print self and DC
   call destroy_self(self_new)

!  ==  Compute green function self -> G(k) 
!  ---------------------------------------------------------------------
   call compute_green(cryst_struc,green,paw_dmft,pawang,1,self,opt_self=1,opt_nonxsum=1)

   call integrate_green(cryst_struc,green,paw_dmft,pawang,prtopt=3,opt_ksloc=3,opt_diff=1) !,opt_nonxsum=1)

   call printocc_green(green,5,paw_dmft,3,chtype="DMFT")
!  call printocc_green(green,9,paw_dmft,3,chtype="DMFT FULL")
   if(paw_dmft%lpsichiortho==1.and.paw_dmft%dmft_prgn==1) then
     call local_ks_green(green,paw_dmft,prtopt=1)
   end if

!  ==  Find fermi level
!  ---------------------------------------------------------------------
   call fermi_green(cryst_struc,green,paw_dmft,pawang,self)
!  call leave_new('COLL')

!  ==  Compute Energy with Mixed self-energy and green function  recomputed with new self
!  ---------------------------------------------------------------------
!  green= lattice green function computed from self for a given chemical potential mu (self_mixed,mu)
!  green= local green function is computed from lattice green function(self_mixed,mu)
!  green= occupations are computed with lattice green   function(self_mixed,mu)
   call compute_energy(cryst_struc,energies_dmft,green,paw_dmft,pawprtvol,pawtab,self,occ_type="nlda",part=part3)

!  == Save self on disk
!  ---------------------------------------------------------------------
   call timab(627,1,tsec)
   call rw_self(self,paw_dmft,prtopt=2,opt_rw=2)
   call timab(627,2,tsec)

!  == Test convergency
!  ---------------------------------------------------------------------
   char_enddmft="DMFT (end of DMFT loop)"
   if(green%ifermie_cv==1.and.self%iself_cv==1.and.green%ichargeloc_cv==1.and.paw_dmft%idmftloop>1) then
     write(message,'(a,8x,a)') ch10,"DMFT Loop is converged !"
     call wrtout(std_out,message,'COLL')
     char_enddmft="converged DMFT"
     exit
   end if
!  =======================================================================
!  === end dmft loop  ====================================================
 end do
!=========================================================================

!== Save self on disk
!-------------------------------------------------------------------------
 call timab(627,1,tsec)
 call rw_self(self,paw_dmft,prtopt=2,opt_rw=2)
 call timab(627,2,tsec)

 !paw_dmft%idmftloop=0

 write(message,'(2a,13x,a)') ch10,' =====  DMFT Loop :  END          ',&
& '========'
 call wrtout(std_out,message,'COLL')

 ! compute Edc for U=1 and J=U/J
 call init_energy(cryst_struc,energies_tmp)
 !call compute_ldau_energy(cryst_struc,energies_tmp,green,paw_dmft,pawtab)
 call compute_ldau_energy(cryst_struc,energies_tmp,green,paw_dmft,pawtab,paw_dmft%forentropyDMFT%J_over_U)
 call data4entropyDMFT_setDc(paw_dmft%forentropyDMFT,energies_tmp%e_dc(:))
 call destroy_energy(energies_tmp,paw_dmft)

!== Compute final values for green functions, occupations, and spectral function
!--------------------------------------------------------------------------------
!Do not compute here, because, one want a energy computed after the
!solver (for Hubbard I and LDA+U).
 call compute_green(cryst_struc,green,paw_dmft,pawang,1,self,opt_self=1,opt_nonxsum=1)
 call integrate_green(cryst_struc,green,paw_dmft,pawang,prtopt=2,opt_ksloc=3) !,opt_nonxsum=1)
!call compute_energy(cryst_struc,energies_dmft,green,paw_dmft,pawprtvol,pawtab,self,opt=0)
 idmftloop=paw_dmft%idmftloop
 paw_dmft%idmftloop=0
 call compute_energy(cryst_struc,energies_dmft,green,paw_dmft,pawprtvol,pawtab,self,occ_type="nlda",part="band")
 paw_dmft%idmftloop=idmftloop

!write(message,'(2a,13x,a)') ch10,' =====  DMFT Loop is finished'
!call wrtout(ab_out,message,'COLL')
!write(std_out,*) "PRINTOCC INITIAL"
 call printocc_green(green,9,paw_dmft,3,chtype=char_enddmft)
!write(std_out,*) "KS=czero"
!green%occup%ks=czero
!write(std_out,*) "PRINTOCC AFTER KS=0"
!call printocc_green(green,9,paw_dmft,3,chtype="converged DMFT")
!write(std_out,*) "UPFOLD_OPER"
!call upfold_oper(green%occup,paw_dmft,1)
!write(std_out,*) "PRINTOCC AFTER UPFOLD_OPER"
!call printocc_green(green,9,paw_dmft,3,chtype="converged DMFT")
!write(std_out,*) "MATLU=czero"
!green%occup%matlu(1)%mat=czero
!green%occup%ks(:,:,:,:)=cmplx(real(green%occup%ks(:,:,:,:)))
!write(std_out,*) "PRINTOCC AFTER MATLU=0 AND IMAG KS=0"
!call printocc_green(green,9,paw_dmft,3,chtype="converged DMFT")
!write(std_out,*) "LOC_OPER"
!call loc_oper(green%occup,paw_dmft,1)
!write(std_out,*) "PRINTOCC AFTER LOC_OPER"
!call printocc_green(green,9,paw_dmft,3,chtype="converged DMFT")
!call flush_unit(std_out)
!call leave_new('COLL')
 if(paw_dmft%dmft_solv<=2.and.paw_dmft%prtdos>=1) then
   call spectral_function(cryst_struc,green,hu,&
&   paw_dmft,pawang,pawtab,self,pawprtvol) 
 end if
 call destroy_green(weiss)
 call destroy_green(green)
!todo_ab rotate back density matrix into unnormalized basis just for
!printout 
 call destroy_hu(hu,cryst_struc%ntypat,paw_dmft%dmftqmc_t2g)
 call destroy_self(self)
 call destroy_energy(energies_dmft,paw_dmft)

 write(message,'(2a,13x,a)') ch10,' =====  DMFT  :  END          ',&
& '========'
 call wrtout(std_out,message,'COLL')

 ABI_DATATYPE_DEALLOCATE(hu)

 DBG_EXIT("COLL")

end subroutine dmft_solve
!!***

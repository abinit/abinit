!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dmft
!! NAME
!!  m_dmft
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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

MODULE m_dmft

 use defs_basis
 implicit none

 private

 public :: dmft_solve
 public :: impurity_solve
 public :: dyson
 public :: spectral_function
!!***

contains

!!****f* ABINIT/dmft_solve
!! NAME
!! dmft_solve
!!
!! FUNCTION
!! Solve the DMFT loop from PAW data.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2019 ABINIT group (BAmadon)
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
!!      m_spectral_function,timab,wrtout
!!
!! SOURCE

subroutine dmft_solve(cryst_struc,istep,lda_occup,paw_dmft,pawang,pawtab,pawprtvol)


 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_abicore
 use m_data4entropyDMFT

 use m_time,           only : timab
 use m_pawang, only : pawang_type
 use m_pawtab, only : pawtab_type
 use m_paw_dmft, only: paw_dmft_type
 use m_crystal, only : crystal_t
 use m_green, only : green_type, destroy_green, icip_green,init_green,&
&                    print_green,printocc_green,&
&                    integrate_green,copy_green,compute_green,&
&                    check_fourier_green,local_ks_green,fermi_green
 use m_oper, only : oper_type,diff_oper,upfold_oper,loc_oper!,upfold_oper,init_oper,destroy_oper,print_oper
 use m_self, only : self_type,initialize_self,destroy_self,print_self,dc_self,rw_self,new_self,make_qmcshift_self
 use m_hu, only : hu_type,init_hu,destroy_hu
 use m_energy, only : energy_type,init_energy,destroy_energy,compute_energy,print_energy,compute_ldau_energy
 use m_matlu, only : print_matlu,sym_matlu!,identity_matlu
 use m_datafordmft, only : psichi_renormalization
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
 integer :: itypat,natomcor,iatom
 logical :: etot_var
 character(len=200) :: char_enddmft
! type
 type(green_type) :: green
 type(green_type) :: greenlda
 type(hu_type),allocatable :: hu(:)
 type(green_type) :: weiss
 type(self_type) :: self
 type(self_type) :: self_new
 !type(oper_type) :: self_minus_hdc_oper
 type(energy_type) :: energies_dmft
 type(energy_type) :: energies_tmp
 character(len=500) :: message
 character(len=5) :: thdyn
 character(len=4) :: part2,part3
!************************************************************************

 DBG_ENTER('COLL')
 my_rank = xmpi_comm_rank(paw_dmft%spacecomm)

 check=paw_dmft%dmftcheck ! checks enabled
 !paw_dmft%dmft_fermi_prec=tol5
 paw_dmft%dmft_fermi_prec=paw_dmft%dmft_charge_prec*ten
!paw_dmft%dmft_charge_prec=20_dp ! total number of electron.
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
 natomcor=0
 do iatom=1,paw_dmft%natom
   if(paw_dmft%lpawu(iatom).ne.-1) then
     natomcor=natomcor+1
   end if
 end do
 opt_renorm=3
! write(6,*) "natomcor",natomcor
! if(natomcor>1) opt_renorm=2
 if(paw_dmft%nspinor==2.and.paw_dmft%dmft_solv==8) opt_renorm=2 ! necessary to use hybri_limit in qmc_prep_ctqmc
                                                                ! ought to be  generalized  in the  future
 if(paw_dmft%dmft_solv/=-1) then
   call psichi_renormalization(cryst_struc,paw_dmft,pawang,opt=opt_renorm)


   ! check that loc_oper(upfold_oper)=I
   !  call init_oper(paw_dmft,self_minus_hdc_oper)
   !  call identity_matlu(self_minus_hdc_oper%matlu,paw_dmft%natom)
   !  call print_matlu(self_minus_hdc_oper%matlu,paw_dmft%natom,1)
   !  call upfold_oper(self_minus_hdc_oper,paw_dmft,1)
   !  call print_oper(self_minus_hdc_oper,9,paw_dmft,3)
   !  call loc_oper(self_minus_hdc_oper,paw_dmft,1)
   !  call sym_matlu(cryst_struc,self_minus_hdc_oper%matlu,pawang)
   !  call print_matlu(self_minus_hdc_oper%matlu,paw_dmft%natom,1)
   !  call destroy_oper(self_minus_hdc_oper)


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
 call init_hu(cryst_struc,pawtab,hu,paw_dmft%dmftqmc_t2g,paw_dmft%dmftqmc_x2my2d)
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
!   call abi_abort('COLL')

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
!   call abi_abort('COLL')

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
!  call abi_abort('COLL')

!  ==  Compute Energy with Mixed self-energy and green function  recomputed with new self
!  ---------------------------------------------------------------------
!  green= lattice green function computed from self for a given chemical potential mu (self_mixed,mu)
!  green= local green function is computed from lattice green function(self_mixed,mu)
!  green= occupations are computed with lattice green   function(self_mixed,mu)
   call compute_energy(cryst_struc,energies_dmft,green,paw_dmft,pawprtvol,pawtab,self,occ_type="nlda",part=part3)

!  == Save self on disk
!  ---------------------------------------------------------------------
   call timab(627,1,tsec)
   call rw_self(self,paw_dmft,prtopt=2,opt_rw=2,pawang=pawang,cryst_struc=cryst_struc)
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
!call abi_abort('COLL')
 if(paw_dmft%dmft_solv<=2.and.paw_dmft%prtdos>=1) then
   call spectral_function(cryst_struc,green,hu,paw_dmft,pawang,pawtab,self,pawprtvol)
 end if
 call destroy_green(weiss)
 call destroy_green(green)
!todo_ab rotate back density matrix into unnormalized basis just for
!printout
 call destroy_hu(hu,cryst_struc%ntypat,paw_dmft%dmftqmc_t2g,paw_dmft%dmftqmc_x2my2d)
 call destroy_self(self)
 call destroy_energy(energies_dmft,paw_dmft)

 write(message,'(2a,13x,a)') ch10,' =====  DMFT  :  END          ',&
& '========'
 call wrtout(std_out,message,'COLL')

 ABI_DATATYPE_DEALLOCATE(hu)

 DBG_EXIT("COLL")

end subroutine dmft_solve
!!***


!!****f* ABINIT/impurity_solve
!! NAME
!! impurity_solve
!!
!! FUNCTION
!! Solve the Impurity problem
!!
!! COPYRIGHT
!! Copyright (C) 1999-2019 ABINIT group (BAmadon)
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
!!      m_dmft
!!
!! CHILDREN
!!      copy_green,destroy_green_tau,flush_unit,fourier_green,m_hubbard_one
!!      init_green_tau,integrate_green,ldau_self,print_green,print_matlu
!!      printocc_green,qmc_prep_ctqmc,timab,trace_oper,wrtout
!!
!! SOURCE

subroutine impurity_solve(cryst_struc,green,hu,paw_dmft,&
& pawang,pawtab,self_old,self_new,weiss,pawprtvol)


 use defs_basis
 use defs_abitypes
 use m_errors
 use m_abicore

 use m_time,    only : timab
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
 use m_hubbard_one, only : hubbard_one
 use m_ldau_self, only : ldau_self
 use m_forctqmc, only : qmc_prep_ctqmc
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


!!****f* ABINIT/dyson
!! NAME
!! dyson
!!
!! FUNCTION
!! Use the Dyson Equation to compute self-energy from green function
!!
!! COPYRIGHT
!! Copyright (C) 1999-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  istep    =  step of iteration for LDA.
!!  lda_occup
!!  mpi_enreg=informations about MPI parallelization
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!  opt_weissgreen= 1: compute weiss from green and self
!!                = 2: compute green from weiss and self
!!                = 4: compute green from weiss and self without
!!                     inversion of weiss  (weiss is previously inverted)
!!
!! OUTPUT
!!  edmft  = energy in DMFT.
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!
!! NOTES
!!
!! PARENTS
!!      dmft_solve
!!
!! CHILDREN
!!      add_matlu,copy_green,destroy_green,init_green,inverse_oper,timab,wrtout
!!
!! SOURCE

subroutine dyson(green,paw_dmft,self,weiss,opt_weissself)

 use defs_basis
 use m_abicore
 use m_errors

 use m_time,         only : timab
 use m_paw_dmft, only: paw_dmft_type
 use m_crystal, only : crystal_t
 use m_green, only : green_type, destroy_green,init_green,copy_green
 use m_oper, only : oper_type,inverse_oper
 use m_matlu, only : matlu_type,add_matlu,print_matlu
 use m_self, only : self_type
 implicit none

!Arguments ------------------------------------
!scalars
 type(green_type),intent(inout)  :: green
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(self_type),intent(inout)  :: self
 type(green_type),intent(inout)  :: weiss
 integer,intent(in) :: opt_weissself
! type(paw_dmft_type), intent(inout)  :: paw_dmft

!Local variables ------------------------------
!arrays
 real(dp) :: tsec(2)
!scalars
 integer :: ifreq,natom,nspinor,nsppol,weissinv
 type(green_type)  :: greeninv
 character(len=500) :: message
! type
! type(matlu_type), pointer :: matlutemp,matlu1,matlu2
!************************************************************************

 call timab(623,1,tsec)
 DBG_ENTER("COLL")
 natom=green%oper(1)%natom
 nsppol=green%oper(1)%nsppol
 nspinor=green%oper(1)%nspinor
 weissinv=1
 if(opt_weissself==1) then
   write(message,'(2a,i3,13x,a)') ch10,'  ===  Use Dyson Equation => weiss '
   call wrtout(std_out,message,'COLL')
 else if(opt_weissself==2) then
   write(message,'(2a,i3,13x,a)') ch10,'  ===  Use Dyson Equation => self'
   call wrtout(std_out,message,'COLL')
 end if

!call xmpi_barrier(spaceComm)
 if(paw_dmft%dmft_solv==2) weissinv=0
 call init_green(greeninv,paw_dmft,opt_oper_ksloc=2,wtype=green%w_type)
 call copy_green(green,greeninv,opt_tw=2)

 do ifreq=1,green%nw
   if(opt_weissself==1) then
     call inverse_oper(greeninv%oper(ifreq),option=1,prtopt=1)
!    warning green is now inversed
     call add_matlu(greeninv%oper(ifreq)%matlu,self%oper(ifreq)%matlu,&
&     weiss%oper(ifreq)%matlu,natom,sign_matlu2=1)
     call inverse_oper(weiss%oper(ifreq),option=1,prtopt=1)
   else if(opt_weissself==2) then

   ! write(59,*) paw_dmft%omega_lo(ifreq), real(green%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(green%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
   ! write(61,*) paw_dmft%omega_lo(ifreq), real(greeninv%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(greeninv%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
   ! write(60,*) paw_dmft%omega_lo(ifreq), real(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
!    call inverse_oper(weiss%oper(ifreq),option=1,prtopt=1)
!    write(62,*) paw_dmft%omega_lo(ifreq), real(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
!    call inverse_oper(weiss%oper(ifreq),option=1,prtopt=1)
!    write(63,*) paw_dmft%omega_lo(ifreq), real(self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))

!    write(std_out,*) "-----------------------IFREQ",ifreq
!    call print_matlu(greeninv%oper(ifreq)%matlu,paw_dmft%natom,1,opt_diag=-1)
     call inverse_oper(greeninv%oper(ifreq),option=1,prtopt=1)
!    call print_matlu(greeninv%oper(ifreq)%matlu,paw_dmft%natom,1,opt_diag=-1)
!    If paw_dmft%dmft_solv==2, then inverse of weiss function is
!    computed in m_hubbard_one.F90
     if(weissinv/=0) then
       call inverse_oper(weiss%oper(ifreq),option=1,prtopt=1)
     end if

!    write(std_out,*) weiss%oper(1)%matlu(ifreq)%mat(1,1,1,1,1),"-",greeninv%oper(ifreq)
     call add_matlu(weiss%oper(ifreq)%matlu,greeninv%oper(ifreq)%matlu,&
&     self%oper(ifreq)%matlu,natom,sign_matlu2=-1)

   ! write(64,*) paw_dmft%omega_lo(ifreq), real(greeninv%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(greeninv%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
   ! write(65,*) paw_dmft%omega_lo(ifreq), real(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
   ! write(66,*) paw_dmft%omega_lo(ifreq), real(self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
   else
     message = " BUG in dyson.F90"
     MSG_BUG(message)
   end if
 end do

 call destroy_green(greeninv)


 call timab(623,2,tsec)
 DBG_EXIT("COLL")

end subroutine dyson
!!***

!!****f* m_dmft/spectral_function
!! NAME
!! spectral_function
!!
!! FUNCTION
!! Print the spectral function computed from Green function in real frequency
!!
!! COPYRIGHT
!! Copyright (C) 1999-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  green  <type(green_type)>= green function data
!!  hu  <type(hu_type)>= datatype of type hu
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  self <type(self_type)>= variables related to self-energy
!!  prtopt= option for printing
!!
!! OUTPUT
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!
!! NOTES
!!
!! PARENTS
!!      m_dmft
!!
!! CHILDREN
!!      compute_green,copy_green,copy_matlu,dc_self,destroy_green,destroy_self
!!      dyson,hubbard_one,init_green,initialize_self,ldau_self,print_green
!!      rw_self,wrtout
!!
!! SOURCE

subroutine spectral_function(cryst_struc,green,hu,paw_dmft,&
& pawang,pawtab,self_old,prtopt)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_abicore

 use m_crystal, only : crystal_t
 use m_green, only : init_green,green_type,print_green,copy_green,compute_green,destroy_green
 use m_matlu, only : copy_matlu
 use m_paw_dmft, only : paw_dmft_type
 use m_hu, only : hu_type
 use m_self, only : self_type,initialize_self,dc_self,destroy_self,rw_self
 use m_energy, only : energy_type
 use m_pawang, only : pawang_type
 use m_pawtab, only : pawtab_type
 use m_hubbard_one, only : hubbard_one
 use m_ldau_self, only : ldau_self
 implicit none

!Arguments ------------------------------------
!scalars
! type(pawang_type), intent(in) :: pawang
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type), intent(in) :: green
 type(hu_type),intent(inout) :: hu(cryst_struc%ntypat)
 !type(MPI_type), intent(inout) :: mpi_enreg
 type(pawang_type), intent(in) :: pawang
 type(pawtab_type),intent(in)  :: pawtab(cryst_struc%ntypat)
 type(self_type), intent(inout) :: self_old
 type(paw_dmft_type), intent(inout)  :: paw_dmft
 integer, intent(in) :: prtopt

!Local variables ------------------------------
 character(len=500) :: message
 type(green_type) :: greenr
 type(green_type) :: weissr
 type(self_type) :: selfr
!scalars
!************************************************************************
!character(len=500) :: message

!   opt_oper_ksloc=3 to be able to compute spectral function.
 call init_green(greenr,paw_dmft,opt_oper_ksloc=3,wtype="real")
 call init_green(weissr,paw_dmft,wtype="real")
 call copy_matlu(green%occup%matlu,greenr%occup%matlu,paw_dmft%natom)
 call initialize_self(selfr,paw_dmft,wtype="real")
!=======================================================================
!== Solve impurity model with green function for real frequency
!=======================================================================
 write(message,'(2a,i3,13x,a)') ch10,'  ===  Write Spectral function'
 call wrtout(std_out,message,'COLL')
 if(abs(paw_dmft%dmft_solv)==1) then

!  == LDA+U for test
!  -------------------
   call ldau_self(cryst_struc,greenr,paw_dmft,&
&   pawtab,selfr,opt_ldau=1,prtopt=prtopt)
 else if(abs(paw_dmft%dmft_solv)==2) then

!  == Hubbard One
!  -------------------
   call hubbard_one(cryst_struc,greenr,hu,paw_dmft,&
&   pawang,prtopt,self_old%hdc,weissr)

 else if(abs(paw_dmft%dmft_solv)==4) then

!  == Nothing
!  -------------------
   message = "spectral_function: This section of code is disabled!"
   MSG_ERROR(message)
   call copy_green(weissr,greenr,opt_tw=1)

 else if(abs(paw_dmft%dmft_solv)>=5) then

!  == Nothing
!  -------------------
   MSG_ERROR("Stopping before copy_green")
   call copy_green(weissr,greenr,opt_tw=1)

 else if(abs(paw_dmft%dmft_solv)==0) then

!  == Nothing
!  -------------------
!  weiss%occup%has_operks=0 -> only local part is duplicated
   call copy_green(weissr,greenr,opt_tw=2)
 end if

!=======================================================================
!== Integrate green function and printout occupations
!For dmft_solv=-1,0,or 1 , the green function was not computed: it
!cannot be integrated
!=======================================================================
 call dc_self(green%charge_matlu_solver,cryst_struc,hu,selfr,paw_dmft%dmft_dc,prtopt)
 if(abs(paw_dmft%dmft_solv)/=1.and.paw_dmft%dmft_solv/=0) then
   call dyson(greenr,paw_dmft,selfr,weissr,opt_weissself=2)
 end if
 call compute_green(cryst_struc,greenr,paw_dmft,pawang,1,selfr,opt_self=1)
 call print_green("realw",greenr,4,paw_dmft,pawprtvol=3)
 call rw_self(selfr,paw_dmft,prtopt=2,opt_rw=2)

 if(abs(prtopt)>0) then
 end if
 call destroy_self(selfr)
 call destroy_green(weissr)
 call destroy_green(greenr)

end subroutine spectral_function
!!***


END MODULE m_dmft
!!***

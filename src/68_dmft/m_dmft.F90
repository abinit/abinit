!!****m* ABINIT/m_dmft
!! NAME
!!  m_dmft
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2025 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

! nvtx related macro definition
#include "nvtx_macros.h"

MODULE m_dmft

 use defs_abitypes
 use defs_basis
 !use netcdf
 use m_xmpi
 use m_abicore
 use m_data4entropyDMFT
 use m_errors

 use m_crystal, only : crystal_t
 use m_datafordmft, only : chipsi_print,chipsi_renormalization,compute_wannier,print_wannier
 use m_dftu_self, only : dftu_self
 use m_energy, only : compute_dftu_energy,compute_energy,compute_free_energy,&
                    & destroy_energy,energy_type,init_energy
 use m_forctqmc, only : ctqmc_calltriqs_c,qmc_prep_ctqmc
 use m_green, only : check_fourier_green,compute_green,copy_green,destroy_green,destroy_green_tau, &
                   & fermi_green,fourier_green,green_type,icip_green,init_green,init_green_tau,integrate_green, &
                   & local_ks_green,print_green,printocc_green
 use m_hu, only : destroy_hu,hu_type,init_hu
 use m_hubbard_one, only : hubbard_one
 use m_matlu, only : add_matlu,copy_matlu,destroy_matlu,diff_matlu,identity_matlu,init_matlu,inverse_matlu, &
                   & matlu_type,print_matlu,sym_matlu,xmpi_matlu
 use m_oper, only : destroy_oper,diff_oper,downfold_oper,gather_oper,init_oper,inverse_oper,oper_type,trace_oper
 use m_paw_dmft, only : paw_dmft_type
 use m_pawang, only : pawang_type
 use m_pawtab, only : pawtab_type
 use m_self, only : dc_self,destroy_self,initialize_self,new_self,print_self,rw_self,self_type
 use m_time, only : timab

#ifdef HAVE_GPU_MARKERS
 use m_nvtx_data
#endif

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
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  istep = iteration step of the DFT+DMFT self-consistent cycle.
!!  dft_occup <type(oper_type)> = DFT occupations numbers of the correlated orbitals
!!  mpi_enreg=information about MPI parallelization
!!  paw_dmft <type(paw_dmft_type)> =  data for self-consistent DFT+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pawprtvol  = option for printing
!!
!! OUTPUT
!!  paw_dmft <type(paw_dmft_type)> = data for self-consistent DFT+DMFT calculations.
!!
!! NOTES
!!
!! SOURCE

subroutine dmft_solve(cryst_struc,istep,dft_occup,mpi_enreg,paw_dmft,pawang,pawtab,pawprtvol)

!Arguments ------------------------------------
 integer, intent(in) :: istep,pawprtvol
 type(MPI_type), intent(in) :: mpi_enreg
 type(pawang_type), intent(in) :: pawang
 type(crystal_t), intent(in) :: cryst_struc
 type(paw_dmft_type), intent(inout)  :: paw_dmft
 type(pawtab_type), intent(inout) :: pawtab(paw_dmft%ntypat)
 type(oper_type), intent(in) :: dft_occup
!Local variables ------------------------------
 integer :: check,dmft_iter,idmftloop,istep_iter,itypat,myproc,natom
 integer :: ntypat,opt_diff,opt_maxent,opt_moments,opt_renorm,prtopt
 !logical :: etot_var
 logical :: dmft_optim,t2g,x2my2d
 real(dp) :: tsec(2)
 character(len=200) :: char_enddmft
 type(green_type) :: green,greendft,weiss
 type(self_type) :: self,self_new
 type(energy_type) :: energies_dmft,energies_tmp
 type(oper_type) :: identity_oper,oper_tmp
 type(hu_type), allocatable :: hu(:)
 character(len=4) :: part2,part3
 !character(len=5) :: thdyn
 character(len=500) :: message
!************************************************************************

 DBG_ENTER('COLL')
 ABI_NVTX_START_RANGE(NVTX_DMFT_SOLVE)

 myproc = paw_dmft%myproc
 check  = paw_dmft%dmftcheck ! checks enabled
 t2g    = (paw_dmft%dmft_t2g == 1)
 x2my2d = (paw_dmft%dmft_x2my2d == 1)
 natom  = paw_dmft%natom
 ntypat = paw_dmft%ntypat
 dmft_iter  = paw_dmft%dmft_iter
 opt_maxent = paw_dmft%dmft_prt_maxent
 dmft_optim  = (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7) ! skip some unneeded calls to compute_green with TRIQS
 !paw_dmft%dmft_fermi_prec=tol5
 !paw_dmft%dmft_fermi_prec = paw_dmft%dmft_charge_prec * ten
!paw_dmft%dmft_charge_prec=20_dp ! total number of electron.
 !paw_dmft%dmft_prgn=1
 paw_dmft%dmft_prgn = 0
 !etot_var = .true.
 !thdyn="fcalc"
 !thdyn = "ecalc"
 !if (thdyn == "ecalc") then ! valid
 part2 = "both"
 part3 = "none"
 !else if (thdyn == "fcalc") then ! not tested
 !  part2 = "corr"
 !  part3 = "band"
 !end if

 opt_moments = merge(1,0,paw_dmft%dmft_solv==6.or.paw_dmft%dmft_solv==7)
 prtopt = merge(2,0,dmft_optim)
 opt_diff = merge(1,0,dmft_optim)

 if (check == 1) then
   write(message,'(2a)') ch10,' DMFT Checks are enabled '
 else
   write(message,'(2a)') ch10,' DMFT Checks will not be performed'
 end if ! check
 call wrtout(std_out,message,'COLL')

 if (istep == 0) then
   message = ' istep should not be equal to zero'
   ABI_BUG(message)
 end if

 !spaceComm=paw_dmft%spacecomm
 !if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_kpt
 !call xmpi_barrier(spaceComm)

 call initialize_self(self,paw_dmft,opt_moments=opt_moments)
 call init_energy(energies_dmft,natom)

!===========================================================================
!==  First construct DFT green function (Init, Compute, Integrate, Print)
!===========================================================================
 write(message,'(6a)') ch10," ==========================================================================", &
                     & ch10," =====  Check: DFT Green's Function Calculation with unnormalized orbitals",&
                     & ch10," =========================================================================="
 call wrtout(std_out,message,'COLL')
 call icip_green("DFT",greendft,paw_dmft,3,self,opt_moments=opt_moments)
 !call print_green('DFT_NOT_renormalized',greendft,1,paw_dmft,pawprtvol=1,opt_wt=1)

!== Compare greendft%occup and dft_occup: check that DFT green function is fine
!----------------------------------------------------------------------
 write(message,'(2a)') ch10," == Compare local occupations from DFT Green's function &
                        &with the downfold of the Fermi-Dirac occupations =="
 call wrtout(std_out,message,'COLL')
 if(paw_dmft%dmft_magnfield .gt. 0) then
   write(message, '(2a,a)') ch10, 'Warning: Check in local occupation is removed due to applied magnetic field'
   call wrtout(std_out,message,'COLL')
 else
   call diff_oper("occupations from DFT Green's function","Fermi-Dirac occupations", &
              & greendft%occup,dft_occup,1,paw_dmft%dmft_tolfreq)
 endif
! write(message,'(2a)') ch10,&
!& '  ***** => Warning : diff_oper is suppressed for test'
! call wrtout(std_out,message,'COLL')
 write(message,'(2a)') ch10,"  ***** => Calculation of DFT Green's function is thus correct ****"
 call wrtout(std_out,message,'COLL')
 call destroy_green(greendft)

!== Orthonormalize chipsi
!----------------------------------------------------------------------
 call timab(621,1,tsec(:))
 !natomcor=0
 !do iatom=1,paw_dmft%natom
 !  if(paw_dmft%lpawu(iatom).ne.-1) then
 !    natomcor=natomcor+1
 !  end if
 !end do
 opt_renorm = merge(2,paw_dmft%dmft_wanorthnorm,paw_dmft%nspinor==2.and.(paw_dmft%dmft_solv==8.or.paw_dmft%dmft_solv==9))

 if (paw_dmft%dmft_solv /= -1) then
   call chipsi_renormalization(paw_dmft,opt=opt_renorm)
   if (paw_dmft%dmft_prtwan == 1) then
     call compute_wannier(paw_dmft,mpi_enreg)
     if (myproc == 0) then
       call print_wannier(paw_dmft,istep)
     end if
     ABI_FREE(paw_dmft%wannier)
   end if ! dmft_prtwan

   write(message,'(2a)') ch10,'  == Check downfold(upfold)=identity =='
   call wrtout(std_out,message,'COLL')

   ! Check that downfold_oper(upfold_oper)=I
   call init_oper(paw_dmft,identity_oper,opt_ksloc=2)
   call init_oper(paw_dmft,oper_tmp,opt_ksloc=2)
   call identity_matlu(identity_oper%matlu(:),natom)
   call downfold_oper(oper_tmp,paw_dmft,procb=paw_dmft%distrib%procb(:),iproc=paw_dmft%distrib%me_kpt,option=4)
   call xmpi_matlu(oper_tmp%matlu(:),natom,paw_dmft%distrib%comm_kpt)
   call sym_matlu(oper_tmp%matlu(:),paw_dmft)
   call diff_matlu("Downfold(Upfold)","Identity",oper_tmp%matlu(:),identity_oper%matlu(:),natom,0,tol4)
   call destroy_oper(oper_tmp)
   call destroy_oper(identity_oper)

!  ===========================================================================
!  ==  re-construct DFT green function with new chipsis
!  ===========================================================================
   write(message,'(6a)') &
    & ch10," ========================================================================", &
    & ch10," =====  Check: DFT Green's Function Calculation with normalized orbitals", &
    & ch10," ========================================================================"
 end if ! dmft_solv/=1
 call timab(621,2,tsec(:))
 call wrtout(std_out,message,'COLL')

 call icip_green("DFT renormalized",greendft,paw_dmft,pawprtvol,self,opt_moments=opt_moments,opt_log=paw_dmft%dmft_triqs_entropy)
 !call print_green('DFT_renormalized',greendft,1,paw_dmft,pawprtvol=1,opt_wt=1)

!== Define Interaction from input upawu and jpawu
!----------------------------------------------------------------------
 ABI_MALLOC(hu,(ntypat))
 call init_hu(hu(:),paw_dmft,pawtab(:))

 call dc_self(greendft%charge_matlu(:,:),self%hdc%matlu(:),hu(:),paw_dmft,pawtab(:),greendft%occup%matlu(:))

 ! Need to store idmftloop and set it to zero to avoid useless print_energy in ab_out
 idmftloop = paw_dmft%idmftloop
 paw_dmft%idmftloop = 0
 call compute_energy(energies_dmft,greendft,paw_dmft,pawprtvol,pawtab(:),self,occ_type=" lda",part='both')
 if (paw_dmft%dmft_triqs_entropy == 1) then
   call compute_free_energy(energies_dmft,paw_dmft,greendft,"band")
 end if
 paw_dmft%idmftloop = idmftloop

 if ((paw_dmft%dmft_prgn == 1) .and. (paw_dmft%lchipsiortho == 1)) then
   call local_ks_green(greendft,paw_dmft,prtopt=1)
 end if
!call printocc_green(greendft,9,paw_dmft,3,chtype="DFT GREEN PSICHI")

 write(message,'(7a)') &
    & ch10,' =============================', &
    & ch10,' =====  Define self-energy', &
    & ch10,' =============================',ch10
 call wrtout(std_out,message,'COLL')

 ! Set Hu in density representation for calculation of entropy if needed...
 if (paw_dmft%dmft_entropy > 0) then
   do itypat=1,ntypat
     if (hu(itypat)%lpawu == -1) cycle
     call data4entropyDMFT_setHu(paw_dmft%forentropyDMFT,itypat,dble(hu(itypat)%udens(:,:)))
   end do ! itypat
 end if ! dmft_entropy=1

!== define self from scratch or file and double counting
!----------------------------------------------------------------------
!-  Self allocated

!-   Read self or do self=hdc
 !if(paw_dmft%dmft_solv==4) then
!  write(std_out,*) "shift before rw_self",self%qmc_shift(1)
 !  call make_qmcshift_self(cryst_struc,hu,self)
 !end if
 call timab(627,1,tsec(:))
 call rw_self(self,paw_dmft,2,opt_rw=1,istep_iter=1000*istep)
 call timab(627,2,tsec(:))

!== If QMC is used,  and self energy is read for file, then
!== one does NOT shifts the self-energy because it was already shifted when writed,
!==  and thus then weiss will be shifted
!----------------------------------------------------------------------
!if(paw_dmft%dmft_solv==4.and.paw_dmft%dmft_rslf==1) &
!&           call make_qmcshift_self(cryst_struc,hu,self)
!if(paw_dmft%dmft_solv==4.and.paw_dmft%dmft_rslf/=1) &
!&           call make_qmcshift_self(cryst_struc,hu,self,apply=.true.)

 call destroy_green(greendft)  ! destroy DFT green function
 call print_self(self,"print_dc",paw_dmft,2)

!===========================================================================
!==  Construct green function with the self-energy.
!===========================================================================
 write(message,'(6a)') &
   & ch10," ===================================================================", &
   & ch10," =====  Green's Function Calculation with input self-energy ========", &
   & ch10," ==================================================================="
 call wrtout(std_out,message,'COLL')
 if (dmft_optim) then
   call init_green(green,paw_dmft,opt_moments=opt_moments)
 else
   call icip_green("DFT+DMFT",green,paw_dmft,pawprtvol,self,opt_self=1,opt_moments=opt_moments)
   !call print_green('beforefermi_green',green,1,paw_dmft,pawprtvol=1,opt_wt=1)
!   call abi_abort('COLL')
 end if

!== Find fermi level
!---------------------------------------------------------------------
!write(message,'(2a,i3,13x,a)') ch10,'   ===  Compute green function from self-energy'

 call fermi_green(green,paw_dmft,self)
 call compute_green(green,paw_dmft,0,self,opt_self=1,opt_nonxsum=1,opt_restart_moments=1)
 call integrate_green(green,paw_dmft,prtopt)

 if (dmft_optim) then
   call printocc_green(green,5,paw_dmft,3,chtype="DFT+DMFT")
 end if

!== define weiss field only for the local quantities (opt_oper=2)
!----------------------------------------------------------------------
! write(std_out,*) "nkpt  befreo init_greenweiss",ifreq,paw_dmft%nkpt
 call init_green(weiss,paw_dmft,opt_oper_ksloc=2,opt_moments=opt_moments,opt_moments_ksloc=2,opt_occup_ksloc=2)
! do ifreq=1,weiss%nw
!   write(std_out,*) "nkpt from weiss1",ifreq,weiss%oper(ifreq)%nkpt
! enddo

!== Check fourier transforms
!----------------------------------------------------------------------
 if (check == 1) then
   call check_fourier_green(cryst_struc,green,paw_dmft)
 end if

!== If QMC is used,  and self energy is not read for file, then
!== one shifts the self-energy, and thus then weiss will be shifted
!== after dyson, in a coherent way regarding qmc_shift and qmc_xmu.
!----------------------------------------------------------------------
!if(paw_dmft%dmft_solv==4.and.paw_dmft%dmft_rslf/=1) &
!&           call make_qmcshift_self(cryst_struc,hu,self,apply=.true.)
!if(paw_dmft%dmft_solv==4) write(std_out,*) "shift after make_qmcshift_self",self%qmc_shift(1)

 write(message,'(6a)') &
   & ch10,' ======================================================', &
   & ch10,' =====  DMFT Loop starts here                  ========', &
   & ch10,' ======================================================'
 call wrtout(std_out,message,'COLL')

 ABI_NVTX_START_RANGE(NVTX_DMFT_SOLVE_LOOP)
!=======================================================================
!===  dmft loop  =======================================================
 do idmftloop=1,dmft_iter
   !paw_dmft%idmftloop=idmftloop
   paw_dmft%idmftloop = paw_dmft%idmftloop + 1
!  =======================================================================
   istep_iter = 1000*istep + idmftloop

   write(message,'(2a,i3,13x,a)') ch10,&
     & ' =====  DMFT Loop : ITER number',paw_dmft%idmftloop,'========'
   call wrtout(std_out,message,'COLL')

!  == Dyson Equation G,self -> weiss(w)
!  ---------------------------------------------------------------------
   call dyson(green,paw_dmft,self,weiss,opt_weissself=1)
!   call print_green('afterDyson',green,1,paw_dmft,pawprtvol=1,opt_wt=1)
!   call abi_abort('COLL')

!  == Printout local "occupations" from weiss field  (useless)
   if (abs(pawprtvol) > 3 .and. opt_moments == 0) then
     call integrate_green(weiss,paw_dmft,2,opt_ksloc=2)
     call printocc_green(weiss,5,paw_dmft,3,opt_weissgreen=1)
   end if

!  ===  Prepare data, solve Impurity problem: weiss(w) -> G(w)
!  ---------------------------------------------------------------------
   call initialize_self(self_new,paw_dmft,opt_moments=opt_moments)

   call impurity_solve(cryst_struc,green,hu(:),paw_dmft,pawang,pawtab(:),self,self_new,weiss,pawprtvol) ! weiss-> green, or self if dmft_solv=1
!  if(paw_dmft%dmft_solv==4)  write(std_out,*) "shift after impurity",self%qmc_shift(1)

!  ==  Compute double counting from charge from green_solver
!  ---------------------------------------------------------------------
   if (green%has_charge_matlu_solver /= 2) green%charge_matlu_solver(:,:) = green%charge_matlu(:,:)

   if (paw_dmft%dmft_solv >= 5) then
     call dc_self(green%charge_matlu_solver(:,:),self_new%hdc%matlu(:),hu(:),paw_dmft,pawtab(:),green%occup_tau%matlu(:))
   else
     call dc_self(green%charge_matlu_solver(:,:),self_new%hdc%matlu(:),hu(:),paw_dmft,pawtab(:),green%occup%matlu(:))
   end if

   if (abs(paw_dmft%dmft_solv) >= 5) then
     call destroy_green_tau(green)
   end if

!  ==  Solve dyson equation. G_imp(w), weiss_imp(w) -> Self_imp(w)
!  ---------------------------------------------------------------------
!  if dmft_solv==1, self is computed previously
   if (abs(paw_dmft%dmft_solv) /= 1) then
     call dyson(green,paw_dmft,self_new,weiss,opt_weissself=2)
   end if
!  do ifreq=1,green%nw
!  call sym_matlu(cryst_struc,self%oper(ifreq)%matlu,pawang)
!  enddo

!  ==  Possibility if imposing self (opt_rw==3)
!  ---------------------------------------------------------------------
   call timab(627,1,tsec(:))
   call rw_self(self_new,paw_dmft,prtopt=2,opt_rw=3,istep_iter=istep_iter)
   call timab(627,2,tsec(:))

!  Print dc computed just before and self computed before in dyson or
!  impurity_solve
   if (abs(pawprtvol) >= 3) then
     write(message,'(2a)') ch10,"  == Old self (before impurity solver)"
     call wrtout(std_out,message,'COLL')
     call print_self(self,"print_dc",paw_dmft,2)
     write(message,'(2a)') ch10,"  == New self (from impurity solver)"
     call wrtout(std_out,message,'COLL')
     call print_self(self_new,"print_dc",paw_dmft,2)
   end if ! abs(pawprtvol)>=3

!  if(paw_dmft%dmft_solv==4) write(std_out,*) "shift before computeenergy ",self%qmc_shift(1)
!  ==  Compute Energy with NEW self-energy and edc from green_solver,
!  new local green function and old occupations for eband
!  fermi level not optimized for this self_energy.
!  ---------------------------------------------------------------------
!  green= local green function and local charge comes directly from solver
!  green= ks green function and occupations comes from old_self
   call compute_energy(energies_dmft,green,paw_dmft,pawprtvol,pawtab(:),self_new,occ_type="nlda",part=part2)
   if (paw_dmft%dmft_triqs_entropy == 1) then
     call compute_free_energy(energies_dmft,paw_dmft,green,"impu",self_new)
   end if

!  ==  Mix new and old self_energies and double countings
!  ---------------------------------------------------------------------
   write(message,'(3a)') ch10,"  == Linear mixing of old and new self-energy and double counting",ch10
   call wrtout(std_out,message,'COLL')
   call new_self(self,self_new,paw_dmft) ! self,self_new => self
   write(message,'(2a)') ch10,"  == After mixing,"
     !print *, " my_rank newself", my_rank,self%oper(1)%matlu(1)%mat(1,1,1,1,1)
   call wrtout(std_out,message,'COLL')
   call print_self(self,"print_dc",paw_dmft,2) ! print self and DC
   call destroy_self(self_new)

!  ==  Compute green function self -> G(k)
!  ---------------------------------------------------------------------
   if (.not. dmft_optim) then
     call compute_green(green,paw_dmft,1,self,opt_self=1,opt_nonxsum=1)
     call integrate_green(green,paw_dmft,3,opt_diff=1) !,opt_nonxsum=1)

     call printocc_green(green,5,paw_dmft,3,chtype="DFT+DMFT")
  !  call printocc_green(green,9,paw_dmft,3,chtype="DMFT FULL")
     if(paw_dmft%lchipsiortho == 1 .and. paw_dmft%dmft_prgn == 1) then
       call local_ks_green(green,paw_dmft,prtopt=1)
     end if
   end if ! dmft_optim=0

!  ==  Find fermi level
!  ---------------------------------------------------------------------
   call fermi_green(green,paw_dmft,self)
   call compute_green(green,paw_dmft,0,self,opt_self=1,opt_nonxsum=1,opt_log=paw_dmft%dmft_triqs_entropy,opt_restart_moments=1)
   call integrate_green(green,paw_dmft,prtopt,opt_diff=opt_diff,opt_ksloc=3,opt_fill_occnd=1)

   if (dmft_optim) then
     call printocc_green(green,5,paw_dmft,3,chtype="DFT+DMFT")
   end if

!  call abi_abort('COLL')

!  ==  Compute Energy with Mixed self-energy and green function recomputed with new self
!  ---------------------------------------------------------------------
!  green= lattice green function computed from self for a given chemical potential mu (self_mixed,mu)
!  green= local green function is computed from lattice green function(self_mixed,mu)
!  green= occupations are computed with lattice green   function(self_mixed,mu)
   call compute_energy(energies_dmft,green,paw_dmft,pawprtvol,pawtab(:),self,occ_type="nlda",part=part3)

!  == Save self on disk
!  ---------------------------------------------------------------------
   call timab(627,1,tsec(:))
   call rw_self(self,paw_dmft,prtopt=2,opt_rw=2,opt_maxent=opt_maxent)
   call timab(627,2,tsec(:))

!  == Test convergency
!  ---------------------------------------------------------------------
   char_enddmft = "DFT+DMFT (end of DMFT loop)"
   if (green%ifermie_cv == 1 .and. self%iself_cv == 1 .and. green%ichargeloc_cv == 1 .and. paw_dmft%idmftloop > 1) then
     write(message,'(a,8x,a)') ch10,"DMFT Loop is converged !"
     call wrtout(std_out,message,'COLL')
     char_enddmft = "converged DMFT"
     exit
   end if
!  =======================================================================
!  === end dmft loop  ====================================================
 end do ! idmftloop
 ABI_NVTX_END_RANGE()
!=========================================================================

!== Save self on disk
!-------------------------------------------------------------------------
 if (.not. dmft_optim) then
   call timab(627,1,tsec(:))
   call rw_self(self,paw_dmft,prtopt=2,opt_rw=2)
   call timab(627,2,tsec(:))
 end if

 !paw_dmft%idmftloop=0

 write(message,'(2a,13x,a)') ch10,' =====  DMFT Loop :  END          ','========'
 call wrtout(std_out,message,'COLL')

 if (paw_dmft%dmft_entropy >= 1) then
   ! compute Edc for U=1 and J=U/J
   call init_energy(energies_tmp,natom)
   !call compute_dftu_energy(cryst_struc,energies_tmp,green,paw_dmft,pawtab)
   call compute_dftu_energy(energies_tmp,green,paw_dmft,pawtab(:),paw_dmft%forentropyDMFT%J_over_U)
   call data4entropyDMFT_setDc(paw_dmft%forentropyDMFT,energies_tmp%e_dc(:))
   call destroy_energy(energies_tmp,paw_dmft)
 end if ! dmft_entropy=1

!== Compute final values for green functions, occupations, and spectral function
!--------------------------------------------------------------------------------
!Do not compute here, because, one want a energy computed after the
!solver (for Hubbard I and DFT+U).
 if (.not. dmft_optim) then
   call compute_green(green,paw_dmft,1,self,opt_self=1,opt_nonxsum=1)
   call integrate_green(green,paw_dmft,2,opt_fill_occnd=1) !,opt_nonxsum=1)
 end if
!call compute_energy(cryst_struc,energies_dmft,green,paw_dmft,pawprtvol,pawtab,self,opt=0)
 idmftloop = paw_dmft%idmftloop
 paw_dmft%idmftloop = 0
 call compute_energy(energies_dmft,green,paw_dmft,pawprtvol,pawtab(:),self,occ_type="nlda",part="band")
 paw_dmft%idmftloop = idmftloop
 if (paw_dmft%dmft_triqs_entropy == 1) then
   call compute_free_energy(energies_dmft,paw_dmft,green,"main",self)
 end if

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
 if (paw_dmft%dmft_solv <= 2 .and. paw_dmft%prtdos >= 1) then
   call spectral_function(cryst_struc,green,hu(:),paw_dmft,pawtab(:),self,pawprtvol)
 end if
 call destroy_green(weiss)
 call destroy_green(green)
!todo_ab rotate back density matrix into unnormalized basis just for
!printout
 call destroy_hu(hu(:),ntypat)
 call destroy_self(self)
 call destroy_energy(energies_dmft,paw_dmft)

 write(message,'(2a,13x,a)') ch10,' =====  DMFT  :  END          ','========'
 call wrtout(std_out,message,'COLL')

 ABI_FREE(hu)

 ABI_NVTX_END_RANGE()
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
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  green  <type(green_type)>= green function data
!!  hu <type(hu_type)>= U interaction
!!  paw_dmft =  data for self-consistent DFT+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  self_old,self_new <type(self_type)>= variables related to self-energy
!!  weiss  <type(green_type)>= weiss function data
!!  pawprtvol = option for printing
!!
!! OUTPUT
!!  paw_dmft =  data for self-consistent DFT+DMFT calculations.
!!
!! NOTES
!!
!! SOURCE

subroutine impurity_solve(cryst_struc,green,hu,paw_dmft,pawang,pawtab,&
                        & self_old,self_new,weiss,pawprtvol)

!Arguments ------------------------------------
!scalars
! type(pawang_type), intent(in) :: pawang
 type(crystal_t), intent(in) :: cryst_struc
 type(green_type), intent(inout) :: green,weiss
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(hu_type), intent(inout) :: hu(paw_dmft%ntypat)
 !type(MPI_type), intent(in) :: mpi_enreg
 type(pawang_type), intent(in) :: pawang
 type(pawtab_type), intent(in) :: pawtab(paw_dmft%ntypat)
 type(self_type), intent(inout) :: self_new,self_old
 integer, intent(in) :: pawprtvol
!Local variables ------------------------------
 real(dp) :: tsec(2)
 character(len=500) :: message
! integer iatom,il,i_nd,isppol,lpawu,im,Nd,nrat,nsweeptot
! real(dp) :: acc,kx
! real(dp), allocatable :: correl(:,:),g0(:,:),gtmp(:,:)
!scalars
!************************************************************************
!character(len=500) :: message

 call timab(622,1,tsec(:))
 ABI_NVTX_START_RANGE(NVTX_DMFT_IMPURITY_SOLVE)
!=======================================================================
!== Prepare data for Hirsch Fye QMC
!== NB: for CTQMC, Fourier Transformation are done inside the CTQMC code
!=======================================================================
 !if(abs(paw_dmft%dmft_solv)==4) then
!  == Initialize weiss and green functions for fourier transformation
!  -------------------------------------------------------------------
 !  write(message,'(2a,i3,13x,a)') ch10,'   ===  Initialize Weiss field G_0(tau)'
 !  call wrtout(std_out,message,'COLL')
 !  call init_green_tau(weiss,paw_dmft)
 !  call init_green_tau(green,paw_dmft)
!  in init_solver

!  == Print weiss function G_0(tau=0-) before computation (really useless check)
!  ------------------------------------------------------------------------------
  ! if(abs(pawprtvol)>3) then
  !   write(message,'(2a,i3,13x,a)') ch10,'   ===  Check G_0(tau=0-) first'
  !   call wrtout(std_out,message,'COLL')
  !   call printocc_green(weiss,6,paw_dmft,3)
  ! end if

!  == Fourier transform of weiss Field
!  ------------------------------------
!  for fourier of KS green functions
!  call fourier_green(cryst_struc,weiss,mpi_enreg,paw_dmft,pawang,pawtab,1)
  ! write(message,'(2a,i3,13x,a)') ch10,'   ===  Inverse Fourier Transform w->t of Weiss Field'
  ! call wrtout(std_out,message,'COLL')
  ! call fourier_green(cryst_struc,weiss,paw_dmft,pawang,opt_ksloc=2,opt_tw=-1)

!  == Print weiss function G2_0(tau=0-)
!  --------------------------------------
  ! call printocc_green(weiss,6,paw_dmft,3,opt_weissgreen=1)

!  for fourier of KS green functions
!  call fourier_green(cryst_struc,weiss,mpi_enreg,paw_dmft,pawang,pawtab,1)
!  == Print G_0(tau) in files
!  ---------------------------
  ! if(paw_dmft%dmft_prgn==1) then
  !   call print_green('weiss',weiss,1,paw_dmft,pawprtvol=1,opt_wt=2)
  ! end if

 if (abs(paw_dmft%dmft_solv) >= 5) then
!  == Initialize  green functions for imaginary times
!  -------------------------------------------------------------------
   write(message,'(2a)') ch10,"  ===  Initialize Green's function G(tau)"
   call wrtout(std_out,message,'COLL')
   call init_green_tau(green,paw_dmft)

 end if
!=======================================================================
!== End preparation of QMC
!=======================================================================

!=======================================================================
!== Solve impurity model   =============================================
!=======================================================================
 write(message,'(2a)') ch10,'  ===  Solving impurity model'
 call wrtout(std_out,message,'COLL')
 if (abs(paw_dmft%dmft_solv) == 1) then

!  == DFT+U for test -> self
!  -------------------
   call dftu_self(cryst_struc,green,paw_dmft,pawtab(:),self_new,opt_dftu=1,prtopt=pawprtvol)

 else if (abs(paw_dmft%dmft_solv) == 2) then

!  == Hubbard One -> green
!  -------------------
   call hubbard_one(cryst_struc,green,hu(:),paw_dmft,pawprtvol,self_old%hdc,weiss)

 !else if(abs(paw_dmft%dmft_solv)==4) then

!  == QMC
!  -------------------
 !  call copy_green(weiss,green,opt_tw=1)
!  call qmc_prep
 !  message = '  ===  QMC not yet distributed '
 !  ABI_ERROR(message)
!   call qmc_prep(cryst_struc,green,hu,mpi_enreg,paw_dmft,pawang&
!&   ,pawprtvol,self_old%qmc_xmu,weiss)

 else if (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7) then

   call ctqmc_calltriqs_c(paw_dmft,green,self_old,hu(:),weiss,self_new,pawprtvol)

 else if (abs(paw_dmft%dmft_solv) >= 5) then

!  == Nothing
!  -------------------
!   call copy_green(weiss,green,opt_tw=1)
!   call copy_green(weiss,green,opt_tw=2)

   call qmc_prep_ctqmc(cryst_struc,green,self_old,hu(:),paw_dmft,pawang,pawprtvol,weiss)


 else if (abs(paw_dmft%dmft_solv) == 0) then

!  == Nothing
!  -------------------
!  weiss%occup%has_operks=0 -> only local part is duplicated
   call copy_green(weiss,green,opt_tw=2)
 end if ! dmft_solv
!call print_green("invWeiss",cryst_struc,weiss,3,paw_dmft,pawtab,2)

!=======================================================================
!== Treat data from HF QMC
!=======================================================================
 if (abs(paw_dmft%dmft_solv) >= 4) then
!  propagate qmc_shift (useful for compute_energy)
   !if(abs(paw_dmft%dmft_solv)==4) then
   !  self_new%qmc_shift(:)=self_old%qmc_shift(:)
   !  self_new%qmc_xmu(:)=self_old%qmc_xmu(:)
   !end if

!  == Print local occupations from G(tau)
!  ---------------------------------------

!  == Fourier back transform of green function G(tau)->G(iw_n) and
!  == compute occupations from g(tau)
!  -------------------------------------------------------------------
   !if(abs(paw_dmft%dmft_solv)==4) then
   !  write(message,'(2a,i3,13x,a)') ch10,'   ===  Direct Fourier Transform t->w of Green Function'
   !  call wrtout(std_out,message,'COLL')
   !  call fourier_green(cryst_struc,green,paw_dmft,&
!&     pawang,opt_ksloc=2,opt_tw=1)
  !   do ifreq=1,green%nw
  !     xx= green%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)
  !     write(112,*) paw_dmft%omega_lo(ifreq),real(one/xx),aimag(one/xx)
  !     write(113,*) paw_dmft%omega_lo(ifreq),real(xx),aimag(xx)
  !   end do
  !   call flush_unit(112)
  !   call flush_unit(113)
!   if(paw_dmft%dmft_solv==5) stop
  !   if(pawprtvol>=3) then
  !     write(message,'(a,2x,a,f13.5)') ch10,&    ! debug
!&      " == Print green function for small freq after fourier " ! debug
  !     call wrtout(std_out,message,'COLL')    ! debug
  !     call print_matlu(green%oper(1)%matlu,paw_dmft%natom,1)    ! debug
  !   end if

   !  write(message,'(2a,i3,13x,a)') ch10,'   INVERSE FOURIER OF G0 SUPPRESSED'
   !  call wrtout(std_out,message,'COLL')
   !end if
   if (abs(paw_dmft%dmft_solv) == 888) then
!  == Back fourier transform of G_0(tau) for compensation (try or comment or improve FT).
!  -------------------------------------------------------------------
     write(message,'(2a)') ch10,'   ===  Direct Fourier transform t->w of Weiss'
     call wrtout(std_out,message,'COLL')
     call fourier_green(cryst_struc,weiss,paw_dmft,opt_ksloc=2,opt_tw=1)

     if (pawprtvol >= 3) then
       write(message,'(a,2x,a,f13.5)') ch10,&    ! debug
         & " == Print weiss function for small freq after fourier " ! debug
       call wrtout(std_out,message,'COLL')    ! debug
       call print_matlu(weiss%oper(1)%matlu(:),paw_dmft%natom,1)    ! debug
     end if ! pawprtvol>=3
     call destroy_green_tau(weiss)
   end if ! dmft_solv=888

!  == Destroy tau part of green
!  -------------------------------------------------------------------
   call trace_oper(green%occup_tau,green%charge_ks,green%charge_matlu_solver(:,:),2)
   green%has_charge_matlu_solver = 2

 end if ! dmft_solv>=5
!=======================================================================
!== End Treat data for QMC
!=======================================================================

!=======================================================================
!== Integrate green function and printout occupations
!=======================================================================
!For dmft_solv=-1,0,or 1, the green function was not yet computed: it
!cannot be integrated
!=======================================================================
 if (paw_dmft%dmft_solv >= 2 .and. green%w_type == "imag") then
!  ==  Integrate G(iw_n)
!  ---------------------
   write(message,'(2a)') ch10,"   ===  Integrate local part of Green's function"
   call wrtout(std_out,message,'COLL')
   call integrate_green(green,paw_dmft,2,opt_ksloc=2,opt_after_solver=1)

!  == Print local occupations from integration of G(iw_n)
!  --------------------------------------------------------
   call printocc_green(green,5,paw_dmft,3)

!  == Print G_loc(w)
!  --------------------------------------------------------
   if (paw_dmft%dmft_prgn == 1) then
     call print_green('DMFT_IMPURITY',green,1,paw_dmft,opt_wt=1)
   end if
 end if ! dmft_solv>=2 and w_type="imag"
!stop

 !if(abs(pawprtvol)>0) then
 !end if


 ABI_NVTX_END_RANGE()
 call timab(622,2,tsec(:))

end subroutine impurity_solve
!!***

!!****f* ABINIT/dyson
!! NAME
!! dyson
!!
!! FUNCTION
!! Use the Dyson Equation to compute self-energy from green function
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  paw_dmft = data for self-consistent DFT+DMFT calculations.
!!  self <type(self_type)>= variables related to self-energy
!!  weiss  <type(green_type)>= Weiss field
!!  opt_weissself = 1: compute weiss from green and self
!!                = 2: compute self from green and weiss
!!
!! OUTPUT
!!
!! NOTES
!!
!! SOURCE

subroutine dyson(green,paw_dmft,self,weiss,opt_weissself)

!Arguments ------------------------------------
!scalars
 type(green_type), intent(in) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(self_type), intent(inout) :: self
 type(green_type), intent(inout) :: weiss
 integer, intent(in) :: opt_weissself
! type(paw_dmft_type), intent(inout)  :: paw_dmft
!Local variables ------------------------------
 integer :: ifreq,myproc,natom,nspinor,nsppol,weissinv
 logical :: triqs
 real(dp) :: tsec(2)
 type(matlu_type), allocatable :: greeninv(:)
 character(len=500) :: message
! type
! type(matlu_type), pointer :: matlutemp,matlu1,matlu2
!************************************************************************

 call timab(623,1,tsec(:))
 DBG_ENTER("COLL")

 myproc   = paw_dmft%myproc
 natom    = paw_dmft%natom
 nsppol   = paw_dmft%nsppol
 nspinor  = paw_dmft%nspinor
 triqs    = (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7)
 weissinv = merge(0,1,paw_dmft%dmft_solv==2.or.triqs)

 if (opt_weissself == 1) then
   write(message,'(2a)') ch10,"  ===  Use Dyson's Equation => weiss"
   call wrtout(std_out,message,'COLL')
 else if (opt_weissself == 2) then
   write(message,'(2a)') ch10,"  ===  Use Dyson's Equation => self"
   call wrtout(std_out,message,'COLL')
 end if ! opt_weisself

!call xmpi_barrier(spaceComm)

 ABI_MALLOC(greeninv,(natom))
 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),greeninv(:))

 do ifreq=1,green%nw

   if (green%distrib%procf(ifreq) /= myproc) cycle

   call copy_matlu(green%oper(ifreq)%matlu(:),greeninv(:),natom)
   call inverse_matlu(greeninv(:),natom)

   if (opt_weissself == 1) then

!    warning green is now inversed
     call add_matlu(greeninv(:),self%oper(ifreq)%matlu(:),weiss%oper(ifreq)%matlu(:),natom,1)
     if (.not. triqs) then
       call inverse_oper(weiss%oper(ifreq),2)
     end if

   else if (opt_weissself == 2) then

   ! write(59,*) paw_dmft%omega_lo(ifreq), real(green%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(green%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
   ! write(61,*) paw_dmft%omega_lo(ifreq), real(greeninv%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(greeninv%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
   ! write(60,*) paw_dmft%omega_lo(ifreq), real(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
!    call inverse_oper(weiss%oper(ifreq),option=1,prtopt=1)
!    write(62,*) paw_dmft%omega_lo(ifreq), real(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
!    call inverse_oper(weiss%oper(ifreq),option=1,prtopt=1)
!    write(63,*) paw_dmft%omega_lo(ifreq), real(self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))

!    write(std_out,*) "-----------------------IFREQ",ifreq
!    call print_matlu(greeninv%oper(ifreq)%matlu,paw_dmft%natom,1,opt_diag=-1)
     !call inverse_oper(greeninv%oper(ifreq),option=1,prtopt=1)
!    call print_matlu(greeninv%oper(ifreq)%matlu,paw_dmft%natom,1,opt_diag=-1)
!    If paw_dmft%dmft_solv==2, then inverse of weiss function is
!    computed in m_hubbard_one.F90
     if (weissinv /= 0) then
       call inverse_oper(weiss%oper(ifreq),2)
     end if

!    write(std_out,*) weiss%oper(1)%matlu(ifreq)%mat(1,1,1,1,1),"-",greeninv%oper(ifreq)
     call add_matlu(weiss%oper(ifreq)%matlu(:),greeninv(:),self%oper(ifreq)%matlu(:),natom,-1)

   ! write(64,*) paw_dmft%omega_lo(ifreq), real(greeninv%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(greeninv%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
   ! write(65,*) paw_dmft%omega_lo(ifreq), real(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
   ! write(66,*) paw_dmft%omega_lo(ifreq), real(self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
   else
     message = " BUG in dyson.F90"
     ABI_BUG(message)
   end if ! opt_weissself
 end do ! ifreq

 if (opt_weissself == 1) then
   call gather_oper(weiss%oper(:),weiss%distrib,paw_dmft,opt_ksloc=2)
 else if (opt_weissself == 2) then
   call gather_oper(self%oper(:),self%distrib,paw_dmft,opt_ksloc=2)
 end if ! opt_weisself

 call destroy_matlu(greeninv(:),natom)
 ABI_FREE(greeninv)

 call timab(623,2,tsec(:))
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
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  green  <type(green_type)>= green function data
!!  hu  <type(hu_type)>= datatype of type hu
!!  paw_dmft =  data for self-consistent DFT+DMFT calculations.
!!  self <type(self_type)>= variables related to self-energy
!!  prtopt= option for printing
!!
!! OUTPUT
!!  paw_dmft = data for self-consistent DFT+DMFT calculations.
!!
!! NOTES
!!
!! SOURCE

subroutine spectral_function(cryst_struc,green,hu,paw_dmft,&
& pawtab,self_old,prtopt)

 use m_dftu_self, only : dftu_self
 use m_green, only : compute_green,copy_green,destroy_green,init_green,print_green
 use m_hubbard_one, only : hubbard_one
 use m_matlu, only : copy_matlu
 use m_self, only : dc_self,destroy_self,initialize_self,rw_self

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type), intent(in) :: green
 type(hu_type),intent(inout) :: hu(cryst_struc%ntypat)
 !type(MPI_type), intent(inout) :: mpi_enreg
 type(pawtab_type),intent(inout)  :: pawtab(cryst_struc%ntypat)
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

!  == DFT+U for test
!  -------------------
   call dftu_self(cryst_struc,greenr,paw_dmft,&
&   pawtab,selfr,opt_dftu=1,prtopt=prtopt)
 else if(abs(paw_dmft%dmft_solv)==2) then

!  == Hubbard One
!  -------------------
   call hubbard_one(cryst_struc,greenr,hu,paw_dmft,&
&   prtopt,self_old%hdc,weissr)

 else if(abs(paw_dmft%dmft_solv)==4) then

!  == Nothing
!  -------------------
   message = "spectral_function: This section of code is disabled!"
   ABI_ERROR(message)
   call copy_green(weissr,greenr,opt_tw=1)

 else if(abs(paw_dmft%dmft_solv)>=5) then

!  == Nothing
!  -------------------
   ABI_ERROR("Stopping before copy_green")
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
 call dc_self(green%charge_matlu_solver,selfr%hdc%matlu,hu,paw_dmft,pawtab,green%occup%matlu)
 if(abs(paw_dmft%dmft_solv)/=1.and.paw_dmft%dmft_solv/=0) then
   call dyson(greenr,paw_dmft,selfr,weissr,opt_weissself=2)
 end if
 call compute_green(greenr,paw_dmft,1,selfr,opt_self=1)
 call print_green("realw",greenr,4,paw_dmft)
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

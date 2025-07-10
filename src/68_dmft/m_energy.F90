!!****m* ABINIT/m_energy
!! NAME
!!  m_energy
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

MODULE m_energy

 use defs_basis
 use m_abi_linalg, only : abi_xgemm
 use m_abicore
 use m_errors
 use m_green, only : green_type,occup_fd
 use m_matlu, only : add_matlu,copy_matlu,destroy_matlu,init_matlu, &
                   & matlu_type,trace_prod_matlu
 use m_oper, only : oper_type
 use m_paw_correlations, only : pawuenergy
 use m_paw_dmft, only : paw_dmft_type
 use m_pawtab, only : pawtab_type
 use m_self, only : self_type
 use m_xmpi, only : xmpi_sum

 implicit none

 private

 public :: init_energy
 public :: compute_energy
 public :: compute_migdal_energy
 public :: compute_dftu_energy
 public :: destroy_energy
 public :: print_energy
 public :: compute_noninterentropy
 public :: compute_free_energy
 public :: compute_trace_log_loc
 public :: print_free_energy
!!***

!!****t* m_energy/energy_type
!! NAME
!!  energy_type
!!
!! FUNCTION
!!  This structured datatype contains interaction matrices for the correlated subspace
!!
!! SOURCE

 type, public :: energy_type ! for each typat

  real(dp) :: e_dc_tot

  real(dp) :: e_dcdc

  real(dp) :: e_hu_dftu_tot

  real(dp) :: e_hu_mig_tot

  real(dp) :: e_hu_qmc_tot

  real(dp) :: e_hu_tot

  real(dp) :: eband_dft

  real(dp) :: eband_dmft

  real(dp) :: edmft

  real(dp) :: ekin_imp

  real(dp) :: emig_imp

  real(dp) :: emig_loc

  real(dp) :: fband_dft

  real(dp) :: fband_dmft

  real(dp) :: fband_imp

  real(dp) :: fband_weiss

  real(dp) :: fdmft

  real(dp) :: fimp

  real(dp) :: integral

  !real(dp) :: natom

  real(dp) :: sdmft

  real(dp) :: simp

  real(dp), allocatable :: e_dc(:)

  real(dp), allocatable :: e_hu_dftu(:)

  real(dp), allocatable :: e_hu_mig(:)

  real(dp), allocatable :: e_hu_qmc(:)

  real(dp), ABI_CONTIGUOUS pointer :: e_hu(:) => null()

 end type energy_type

!!***

!----------------------------------------------------------------------


CONTAINS  !========================================================================================
!!***

!!****f* m_energy/init_energy
!! NAME
!! init_energy
!!
!! FUNCTION
!!  Allocate variables used in type energy_type.
!!
!! INPUTS
!!
!! OUTPUTS
!!  energies_dmft  = datastructure for dmft energy
!!  natom = number of atoms
!!
!! SOURCE

subroutine init_energy(energies_dmft,natom)

!Arguments ------------------------------------
!type
 type(energy_type), target, intent(inout) :: energies_dmft
 integer, intent(in) :: natom
!Local variables ------------------------------------
!************************************************************************

 ABI_MALLOC(energies_dmft%e_dc,(natom))
 ABI_MALLOC(energies_dmft%e_hu_dftu,(natom))
 ABI_MALLOC(energies_dmft%e_hu_mig,(natom))
 ABI_MALLOC(energies_dmft%e_hu_qmc,(natom))
 energies_dmft%e_hu => energies_dmft%e_hu_mig(:)
 energies_dmft%e_dc(:)       = zero
 energies_dmft%e_hu_dftu(:)  = zero
 energies_dmft%e_hu_mig(:)   = zero
 energies_dmft%e_hu_qmc(:)   = zero
 energies_dmft%e_dc_tot      = zero
 energies_dmft%e_dcdc        = zero
 energies_dmft%e_hu_dftu_tot = zero
 energies_dmft%e_hu_mig_tot  = zero
 energies_dmft%e_hu_qmc_tot  = zero
 energies_dmft%e_hu_tot      = zero
 energies_dmft%eband_dft     = zero
 energies_dmft%eband_dmft    = zero
 energies_dmft%edmft         = zero
 energies_dmft%ekin_imp      = zero
 energies_dmft%emig_imp      = zero
 energies_dmft%emig_loc      = zero
 energies_dmft%fband_dft     = zero
 energies_dmft%fband_dmft    = zero
 energies_dmft%fband_imp     = zero
 energies_dmft%fband_weiss   = zero
 energies_dmft%fdmft         = zero
 energies_dmft%fimp          = zero
 energies_dmft%integral      = zero
 !energies_dmft%natom         = natom
 energies_dmft%sdmft         = zero
 energies_dmft%simp          = zero

end subroutine init_energy
!!***

!!****f* m_energy/destroy_energy
!! NAME
!! destroy_energy
!!
!! FUNCTION
!!  Deallocate energies_dmft
!!
!! INPUTS
!!  energies_dmft  = datastructure for dmft energy
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_energy(energies_dmft,paw_dmft)

!Arguments ------------------------------------
!scalars
 type(energy_type), intent(inout) :: energies_dmft
 type(paw_dmft_type), intent(inout) :: paw_dmft
!Local variables-------------------------------
! *********************************************************************

 paw_dmft%e_dc  = energies_dmft%e_dc_tot
 paw_dmft%e_hu  = energies_dmft%e_hu_tot
 paw_dmft%sdmft = energies_dmft%sdmft
 paw_dmft%simp  = energies_dmft%simp

 energies_dmft%e_hu => null()
 ABI_SFREE(energies_dmft%e_dc)
 ABI_SFREE(energies_dmft%e_hu_dftu)
 ABI_SFREE(energies_dmft%e_hu_mig)
 ABI_SFREE(energies_dmft%e_hu_qmc)

end subroutine destroy_energy
!!***

!!****f* m_energy/print_energy
!! NAME
!! print_energy
!!
!! FUNCTION
!!  Print different components of DMFT contribution to the internal energy.
!!
!! INPUTS
!!  energies_dmft  = datastructure for dmft energy
!!  pawprtvol = flag for print_energy
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  idmftloop = iteration number of the DFT+DMFT loop
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! SOURCE

subroutine print_energy(energies_dmft,pawprtvol,paw_dmft,idmftloop)

!Arguments ------------------------------------
!type
 type(energy_type), intent(in) :: energies_dmft
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, intent(in) :: pawprtvol,idmftloop
!Local variables-------------------------------
 integer :: iatom,lpawu
 character(len=4) :: tag
 character(len=1000) :: message
! *********************************************************************

 if (abs(pawprtvol) >= 3) then
   do iatom=1,paw_dmft%natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     write(tag,'(i4)') iatom
     write(message,'(a,4x,3a)') ch10,"For Correlated Atom ",trim(adjustl(tag)),","
     call wrtout(std_out,message,'COLL')
     write(message,'(26x,a,f12.6)') " E_hu      =",energies_dmft%e_hu(iatom)
     call wrtout(std_out,message,'COLL')
     write(message,'(26x,a,f12.6)') " E_hu_mig  =",energies_dmft%e_hu_mig(iatom)
     call wrtout(std_out,message,'COLL')
     write(message,'(26x,a,f12.6)') " E_hu_qmc  =",energies_dmft%e_hu_qmc(iatom)
     call wrtout(std_out,message,'COLL')
     write(message,'(26x,a,f12.6)') " E_hu_dftu =",energies_dmft%e_hu_dftu(iatom)
     call wrtout(std_out,message,'COLL')
     write(message,'(26x,a,f12.6)') " E_dc      =",energies_dmft%e_dc(iatom)
     call wrtout(std_out,message,'COLL')
   end do ! iatom
 end if ! abs(pawprtvol)>=3
 write(message,'(a,5x,2a,5x,a,9(a,5x,a,2x,f15.11),a,5x,a)') ch10, &
     & "-----------------------------------------------",ch10, &
     & "--- Energy in DMFT (in Ha)  ",ch10, &
     & "--- E_bandlda (1)  (Ha.) = ",energies_dmft%eband_dft,ch10, &
     & "--- E_banddmft(2)  (Ha.) = ",energies_dmft%eband_dmft,ch10, &
     & "--- E_hu      (3)  (Ha.) = ",energies_dmft%e_hu_tot,ch10, &
     & "--- E_hu_mig  (4)  (Ha.) = ",energies_dmft%e_hu_mig_tot,ch10, &
     & "--- E_hu_qmc  (4)  (Ha.) = ",energies_dmft%e_hu_qmc_tot,ch10, &
     & "--- E_hu_dftu (5)  (Ha.) = ",energies_dmft%e_hu_dftu_tot,ch10, &
     & "--- E_dc      (6)  (Ha.) = ",energies_dmft%e_dc_tot,ch10, &
     & "--- edmft=(    3-6)(Ha.) = ",energies_dmft%edmft,ch10, &
     & "---       (2-1+3-6)(Ha.) = ",energies_dmft%eband_dmft-energies_dmft%eband_dft+energies_dmft%edmft,ch10, &
     & "-----------------------------------------------"
 call wrtout(std_out,message,'COLL')
 if (idmftloop >= 1) then
   write(message,'(a,i3,1x,f15.11,a)') " (Edmft",idmftloop,energies_dmft%edmft,")"
   call wrtout(ab_out,message,'COLL')
 end if

end subroutine print_energy
!!***

!!****f* m_energy/compute_energy
!! NAME
!! compute_energy
!!
!! FUNCTION
!!  Compute and print the different contributions for the DMFT energy.
!!
!! INPUTS
!!  energies_dmft  = datastructure for dmft energy
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawprtvol = flag for printing
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  self  <type(self_type)>= self energy function data
!!  occ_type=  character ("lda" or "nlda") for printing.
!!  part = "band" : compute the DFT and DMFT band energies
!!       = "corr" : compute the interaction energy
!!       = "both" : compute band energy and interaction energy
!!       = "none" : do not print
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  energies_dmft = datastructure for dmft energy
!!
!! SOURCE

subroutine compute_energy(energies_dmft,green,paw_dmft,pawprtvol,pawtab,self,occ_type,part)

!Arguments ------------------------------------
!type
 type(energy_type), target, intent(inout) :: energies_dmft
 type(green_type), intent(in) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(pawtab_type), intent(in) :: pawtab(paw_dmft%ntypat)
 type(self_type), intent(in) :: self
 integer, intent(in) :: pawprtvol
 character(len=4), intent(in) :: occ_type,part
! integer :: prtopt
!Local variables-------------------------------
 integer :: iatom,lpawu
 character(len=500) :: message
! *********************************************************************

 if (part == 'both') then
   if (occ_type == " lda") then
     write(message,'(2a)') ch10,"  == Check: Compute DFT energy terms"
   else
     write(message,'(2a)') ch10,"  == Compute DFT+DMFT energy terms"
   end if
   call wrtout(std_out,message,'COLL')
 else if (part == 'band') then
   write(message,'(2a)') ch10,"  == Compute DFT+DMFT energy terms : Band energy terms"
   call wrtout(std_out,message,'COLL')
 else if (part == 'corr') then
   write(message,'(2a)') ch10,"  == Compute DFT+DMFT energy terms : Correlation energy terms only"
   call wrtout(std_out,message,'COLL')
 !else if(part=='none') then
 end if ! part

! Only imaginary frequencies here
 if (green%w_type == "real" .or. self%w_type == "real") then
   message = 'compute_energy not implemented for real frequency'
   ABI_BUG(message)
 end if

 if (part == 'band' .or. part == 'both') then
   call compute_band_energy(energies_dmft,green,paw_dmft,occ_type,ecalc_dft=1)
 end if
! == Compute Band Energy Alternative version: two steps
!                 == Compute Tr[ln G^{-1}] and -Tr[(Self-hdc)G_dmft]
! -----------------------------------------------------------------------
  ! if (part == 'band') then ! ie if thdyn="fcalc" in m_dmft.F90
!     call compute_B3(cryst_struc,energies_dmft,eband2,green,mpi_enreg,paw_dmft,2,pawang,self,occ_type,0)
!     write(message,'(2a,f10.6)') ch10,"Compute Band energy test KS           ",eband2
!     call wrtout(std_out,message,'COLL')
!     call wrtout(ab_out,message,'COLL')
!
!     call compute_B3(cryst_struc,energies_dmft,eband2,green,mpi_enreg,paw_dmft,2,pawang,self,occ_type,1)
!     write(message,'(2a,f10.6)') ch10,"Compute Band energy test Self statique",eband2
!     call wrtout(std_out,message,'COLL')
!     call wrtout(ab_out,message,'COLL')
!
!! == Compute Band Energy (classical)
!! -----------------------------------------------------------------------
!     call compute_band_energy(energies_dmft,green,paw_dmft,occ_type,fcalc_dft=3)
!     write(std_out,*) paw_dmft%fermie_dft,paw_dmft%fermie
!     write(message,'(2a,f10.6)') ch10,"Compute Band energy ref  free lda -ef ",energies_dmft%eband_dft
!     call wrtout(std_out,message,'COLL')
    ! call compute_band_energy(energies_dmft,green,paw_dmft,occ_type,ecalc_dft=1)
!     write(message,'(2a,f10.6)') ch10,"Compute Band energy ref  -ef          ",energies_dmft%eband_dmft
!     call wrtout(std_out,message,'COLL')
!! call wrtout(ab_out,message,'COLL')
!     write(message,'(2a,f10.6)') ch10,"Compute Band energy ref   lda         ",energies_dmft%eband_dft
!     call wrtout(std_out,message,'COLL')
!! if(occ_type=="nlda") eband2=energies_dmft%eband_dmft
  ! else
   !call compute_band_energy(energies_dmft,green,paw_dmft,occ_type,ecalc_dft=1)
   !endif

 !end if

 if (part == 'corr' .or. part == 'both') then

! == Compute Correlation energy from Migdal formula
! -----------------------------------------------------------------------
   if (occ_type /= " lda") then
     call compute_migdal_energy(energies_dmft%e_hu_mig(:),energies_dmft%e_hu_mig_tot,green,paw_dmft,self)
   end if
! write(std_out,*) "MIGDAL",e_hu_migdal_tot,e_hu_migdal

! == Compute DFT+U interaction energy
! -----------------------------------------------------------------------
   call compute_dftu_energy(energies_dmft,green,paw_dmft,pawtab(:))
   if (abs(paw_dmft%dmft_solv) <= 1) then
     energies_dmft%e_hu => energies_dmft%e_hu_dftu(:)
     energies_dmft%e_hu_tot = energies_dmft%e_hu_dftu_tot
     if ((abs(energies_dmft%e_hu_tot-energies_dmft%e_hu_mig_tot) >= tol6) .and. (occ_type /= " lda")) then
       write(message,'(2a,2e18.8,2x,a)') ch10,'   BUG: Migdal energy and DFT+U energy do not coincide',&
         & energies_dmft%e_hu_tot,energies_dmft%e_hu_mig_tot,occ_type
       ABI_ERROR(message)
     end if
   else if (paw_dmft%dmft_solv == 2 .or. ((paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7) &
     & .and. (.not. paw_dmft%dmft_triqs_measure_density_matrix)) .or. paw_dmft%dmft_solv == 9) then
     energies_dmft%e_hu => energies_dmft%e_hu_mig(:)
     energies_dmft%e_hu_tot = energies_dmft%e_hu_mig_tot
     energies_dmft%e_hu_qmc_tot = energies_dmft%e_hu_tot
   else if (paw_dmft%dmft_solv == 5 .or. paw_dmft%dmft_solv == 8 .or. &
      & ((paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7) .and. &
      & paw_dmft%dmft_triqs_measure_density_matrix) .and. occ_type /= " lda") then
     if (paw_dmft%dmft_solv == 8) then
       write(message,'(2a)') ch10,"Warning, energy is recently computed, not checked"
       call wrtout(std_out,message,'COLL')
     end if
     ! == Compute Correlation energy from QMC correlations.
     ! -----------------------------------------------------------------------
     do iatom=1,paw_dmft%natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       energies_dmft%e_hu_qmc(iatom) = green%ecorr_qmc(iatom)
       energies_dmft%e_hu_qmc_tot = energies_dmft%e_hu_qmc_tot + energies_dmft%e_hu_qmc(iatom)
     end do ! iatom
     energies_dmft%e_hu => energies_dmft%e_hu_qmc(:)
     energies_dmft%e_hu_tot = energies_dmft%e_hu_qmc_tot
   end if ! dmft_solv
!   energies_dmft%edmft=energies_dmft%e_hu_mig_tot-energies_dmft%e_dc_tot
   energies_dmft%edmft = energies_dmft%e_hu_tot - energies_dmft%e_dc_tot

 end if ! part

! if(part='corr'.or.part='both') then
 if (part /= 'none') then
   call print_energy(energies_dmft,pawprtvol,paw_dmft,paw_dmft%idmftloop)
 end if
! write(message,'(2a)') ch10," == The DFT+U self-energy is == "
! call wrtout(std_out,message,'COLL')
! call print_oper(self%oper(1),5,paw_dmft,2)
! a voir: energies_dmft%e_hu_tot = energies_dmft%e_hu_dftu_tot

end subroutine compute_energy
!!***

!!****f* m_energy/compute_band_energy
!! NAME
!! compute_band_energy
!!
!! FUNCTION
!!  Compute the DFT and DMFT band energy
!!
!! INPUTS
!!  energies_dmft  = datastructure for dmft energy
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  occ_type=  character ("lda" or "nlda") for printing.
!!  fcalc_dft= if present, compute free energy/grand potential instead of total energy.
!!           = 2/4 use the DMFT Fermi level
!!           = 2/3 compute grand potential
!!           = 1/3 use the DFT Fermi level
!!  ecalc_dft= 2/4 use the DMFT Fermi level
!!           = 2/3 substract mu*N to the DFT band energy
!!           = 1/3 use the DFT Fermi level
!!  ecalc_dmft= if present, substract mu*N to the DMFT band energy
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  energies_dmft = datastructure for dmft energy
!!
!! SOURCE

subroutine compute_band_energy(energies_dmft,green,paw_dmft,occ_type,ecalc_dft,fcalc_dft,ecalc_dmft)

!Arguments ------------------------------------
 type(energy_type), intent(inout) :: energies_dmft
 type(green_type), intent(in) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 character(len=4), intent(in) :: occ_type
 integer, optional, intent(in) :: ecalc_dft,ecalc_dmft,fcalc_dft
! integer :: prtopt
!Local variables-------------------------------
 integer :: band_index,ib,ibc,ikpt,isppol,nband_k,nkpt,nspinor,nsppol
 real(dp) :: beta,eig,fermie_used,occ,totch2,wtk !,totch3
 character(len=500) :: message
! *********************************************************************

 if (occ_type == " lda") then
   write(message,'(2a)') ch10,"  == Compute DFT Band Energy terms"
 else if (present(fcalc_dft)) then
   write(message,'(2a)') ch10,"  == Compute DFT Free Energy terms"
 else
   write(message,'(2a)') ch10,"  == Compute DMFT Band Energy terms"
 end if
 call wrtout(std_out,message,'COLL')
 beta = one / paw_dmft%temp

! == Compute Band Energy
! -----------------------------------------------------------------------
 if (occ_type == " lda") energies_dmft%eband_dft = zero
 if (.not. present(fcalc_dft)) energies_dmft%eband_dmft = zero
 if (present(fcalc_dft)) energies_dmft%fband_dft = zero
 !totch  = zero
 totch2 = zero
 !totch3 = zero

 nkpt    = paw_dmft%nkpt
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol

 band_index = 0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k = paw_dmft%nband(ikpt+(isppol-1)*nkpt)
     wtk = paw_dmft%wtk(ikpt)
     ibc = 0
     do ib=1,nband_k
       if ((.not. paw_dmft%band_in(ib)) .and. (paw_dmft%dmft_solv /= 6 .and. paw_dmft%dmft_solv /= 7)) cycle
       if (paw_dmft%band_in(ib)) ibc = ibc + 1
       eig = paw_dmft%eigen(ib+band_index)
       if (present(fcalc_dft)) then
         if (fcalc_dft == 1 .or. fcalc_dft == 3) fermie_used = paw_dmft%fermie_dft
         if (fcalc_dft == 2 .or. fcalc_dft == 4) fermie_used = paw_dmft%fermie ! only for B3 terms
         energies_dmft%fband_dft = energies_dmft%fband_dft + wtk*merge(-log(one+exp(-beta*(eig-fermie_used))), &
              & (beta*(eig-fermie_used)-log(one+exp(beta*(eig-fermie_used)))),(eig-fermie_used)>=zero)
       else
         if (occ_type == " lda") then ! usual calculation: total non interacting energy
           fermie_used = paw_dmft%fermie_dft
!            write(std_out,*) "isppol,ikpt,ib",isppol,ikpt,ib
!            write(std_out,*) "paw_dmft%eigen_dft",paw_dmft%eigen_dft(isppol,ikpt,ib)
!            write(std_out,*) green%occup%ks(isppol,ikpt,ib,ib)
!            write(std_out,*) occup_fd(paw_dmft%eigen_dft(isppol,ikpt,ib),paw_dmft%fermie,paw_dmft%temp)
           if (present(ecalc_dft)) then
             if (ecalc_dft == 1 .or. ecalc_dft == 3) fermie_used = paw_dmft%fermie_dft
             if (ecalc_dft == 2 .or. ecalc_dft == 4) fermie_used = paw_dmft%fermie ! only for B3 terms
             occ = occup_fd(eig,fermie_used,paw_dmft%temp)
             if (ecalc_dft == 3 .or. ecalc_dft == 2) then
               energies_dmft%eband_dft = energies_dmft%eband_dft - occ*fermie_used*wtk
               totch2 = totch2 + wtk*occ
             end if
           else
             occ = occup_fd(eig,fermie_used,paw_dmft%temp)
           end if ! present(ecalc_dft)
           energies_dmft%eband_dft = energies_dmft%eband_dft + occ*eig*wtk
         end if ! occ_type=" lda"
         if (paw_dmft%band_in(ib)) then
           occ = dble(green%occup%ks(ibc,ibc,ikpt,isppol))
         else
           occ = occup_fd(eig,paw_dmft%fermie,paw_dmft%temp)
         end if
         energies_dmft%eband_dmft = energies_dmft%eband_dmft + occ*eig*wtk
       !totch3 = totch3 + paw_dmft%wtk(ikpt)*green%occup%ks(ib,ib,ikpt,isppol)
         if (present(ecalc_dmft)) energies_dmft%eband_dmft = energies_dmft%eband_dmft - &
             & occ*paw_dmft%fermie*wtk
       end if ! present(fcalc_dft)
       !end if
     end do ! ib
     band_index = band_index + nband_k
   end do ! ikpt
 end do ! isppol

 if (present(fcalc_dft)) then
   energies_dmft%fband_dft = energies_dmft%fband_dft * paw_dmft%temp
   if (nsppol == 1 .and. nspinor == 1) energies_dmft%fband_dft = energies_dmft%fband_dft * two
   if (fcalc_dft == 1 .or. fcalc_dft == 4) energies_dmft%fband_dft = energies_dmft%fband_dft + &
      & fermie_used*paw_dmft%nelectval
 else
   if (occ_type == " lda" .and. nsppol == 1 .and. nspinor == 1) energies_dmft%eband_dft = two * energies_dmft%eband_dft
   if (nsppol == 1 .and. nspinor == 1) energies_dmft%eband_dmft = two * energies_dmft%eband_dmft
 end if ! present(fcalc_dft)
   !if (fcalc_dft == 3 .or. fcalc_dft == 2) write(std_out,*) "compute_band_energy totch",totch
 !end if

 if (present(ecalc_dft)) then
   if (ecalc_dft == 3 .or. ecalc_dft == 2) write(std_out,*) "compute_band_energy totch2",totch2
 end if
! write(std_out,*) "compute_band_energy totch3",totch3

 if (occ_type == " lda") then
   if (abs(energies_dmft%eband_dft-energies_dmft%eband_dmft) > tol5) then
     write(message,'(5x,3a,15x,a,f12.6,a,15x,a,5x,f12.5)') "Warning: ", &
       & "Differences between band energy with Fermi-Dirac occupations",ch10, &
       & "and occupations from DFT Green's function is:",energies_dmft%eband_dft-energies_dmft%eband_dmft,ch10, &
       & "which is larger than",tol5
     call wrtout(std_out,message,'COLL')
     write(message,'(a)') &
       & "   Action: increase the number of frequencies, or reduce the number of high energy DMFT bands"
     call wrtout(std_out,message,'COLL')
   else
     write(message,'(3a,10x,a,f12.6,a,10x,a,5x,f12.5)')  "          ", &
      & "Differences between band energy with Fermi-Dirac occupations",ch10, &
      & "and occupations from DFT Green's function is:",energies_dmft%eband_dft-energies_dmft%eband_dmft,ch10, &
      & "which is smaller than",tol5
     call wrtout(std_out,message,'COLL')
   end if ! tol
 end if ! occ_type=lda

 if (present(fcalc_dft)) then
   if (abs(energies_dmft%fband_dft-green%trace_log) > tol5) then
     write(message,'(5x,3a,15x,a,f12.6,a,15x,a,5x,f12.5)') "Warning: ", &
       & "Differences between free energy with Fermi-Dirac occupations",ch10, &
       & "and occupations from DFT Green's function is:",energies_dmft%fband_dft-green%trace_log,ch10, &
       & "which is larger than",tol5
     call wrtout(std_out,message,'COLL')
     write(message,'(a)') &
       & "   Action: increase the number of frequencies, or reduce the number of high energy DMFT bands"
     call wrtout(std_out,message,'COLL')
   else
     write(message,'(3a,10x,a,f12.6,a,10x,a,5x,f12.5)')  "          ", &
      & "Differences between free energy with Fermi-Dirac occupations",ch10, &
      & "and occupations from DFT Green's function is:",energies_dmft%fband_dft-green%trace_log,ch10, &
      & "which is smaller than",tol5
     call wrtout(std_out,message,'COLL')
   end if ! tol
 end if ! occ_type=lda

end subroutine compute_band_energy
!!***

!!****f* m_energy/compute_migdal_energy
!! NAME
!! compute_migdal_energy
!!
!! FUNCTION
!!  Computes Midgal energy = 1/2 * Tr(Sigma*G)
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  self  <type(self_type)>= self energy function data
!!  iatom = if present, only computes the contribution of this atom
!!
!! OUTPUT
!!  e_hu_migdal(natom)= Migdal energy for each atom.
!!  e_hu_mig_tot= Total Migdal energy.
!!
!! SOURCE

subroutine compute_migdal_energy(e_hu_migdal,e_hu_migdal_tot,green,paw_dmft,self,iatom)

!#ifdef FC_INTEL
!DEC$ NOOPTIMIZE
!#endif

!Arguments ------------------------------------
 type(green_type), intent(in) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 real(dp), intent(out) :: e_hu_migdal_tot
 real(dp), intent(inout) :: e_hu_migdal(paw_dmft%natom)
 type(self_type), target, intent(in) :: self
 integer, optional, intent(in) :: iatom
! integer :: prtopt
!Local variables-------------------------------
 integer :: i,iatom_,ierr,ifreq,j,myproc,natom,nmoments,nspinor,nsppol,nwlo
 real(dp) :: beta,temp
 complex(dpc) :: omega
 complex(dpc), allocatable :: omega_fac(:),trace_moments(:,:),trace(:)
 type(matlu_type), allocatable :: self_nwlo_re(:)
 type(matlu_type), pointer :: matlu_tmp(:) => null()
 character(len=500) :: message
! *********************************************************************

! Only imaginary frequencies here
 if (green%w_type == "real" .or. self%w_type == "real") then
   message = 'compute_migdal_energy not implemented for real frequency'
   ABI_BUG(message)
 end if

! == Compute Correlation energy from Migdal formula
! -----------------------------------------------------------------------
 beta     = one / paw_dmft%temp
 myproc   = paw_dmft%myproc
 natom    = paw_dmft%natom
 nmoments = 0
 nspinor  = paw_dmft%nspinor
 nsppol   = paw_dmft%nsppol
 nwlo     = green%nw
 temp     = paw_dmft%temp

 if (self%has_moments == 1) nmoments = self%nmoments
 iatom_ = 0
 if (present(iatom)) then
   iatom_ = iatom
   if (self%has_moments == 0) ABI_BUG("You should not be here")
 end if

 if (green%nw /= self%nw) then
   message = 'self and green do not contain the same number of frequencies'
   ABI_BUG(message)
 end if
! write(std_out,*) "beta",beta

 ABI_MALLOC(trace_moments,(natom,nmoments))
 ABI_MALLOC(trace,(natom))

 e_hu_migdal(:) = zero
 trace(:) = czero

 if (self%has_moments == 1) then
   trace_moments(:,:) = czero
   do i=1,nmoments
     do j=1,i
       call trace_prod_matlu(self%moments(j)%matlu(:),green%moments(i-j+1)%matlu(:),natom,trace(:),iatom=iatom_)
       trace_moments(:,i) = trace_moments(:,i) + trace(:)
     end do ! j
   end do ! i
 else
   ABI_MALLOC(self_nwlo_re,(natom))
   ABI_MALLOC(matlu_tmp,(natom))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),matlu_tmp(:))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),self_nwlo_re(:))
   call copy_matlu(self%oper(nwlo)%matlu(:),self_nwlo_re(:),natom,opt_re=1)
 end if ! moments

 do ifreq=1,nwlo

   if (self%distrib%procf(ifreq) /= myproc) cycle

   omega = cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)

   if (self%has_moments == 1) then
     matlu_tmp => self%oper(ifreq)%matlu(:)
   else
     call add_matlu(self%oper(ifreq)%matlu(:),self_nwlo_re(:),matlu_tmp(:),natom,-1)
   end if ! moments

   call trace_prod_matlu(matlu_tmp(:),green%oper(ifreq)%matlu(:),natom,trace(:),iatom=iatom_)

   e_hu_migdal(:) = e_hu_migdal(:) + dble(trace(:))*paw_dmft%wgt_wlo(ifreq)*temp*two

 end do ! ifreq

 call xmpi_sum(e_hu_migdal(:),paw_dmft%spacecomm,ierr)

 ABI_MALLOC(omega_fac,(nmoments))

 do i=1,nmoments
   omega_fac(i) = czero
   do ifreq=nwlo,1,-1 ! NEVER change the summation order and DON'T use the intrinsic SUM
     omega_fac(i) = omega_fac(i) + cone / (paw_dmft%omega_lo(ifreq))**i
   end do
   omega_fac(i) = - two * temp * omega_fac(i) / (j_dpc)**i
   if (i == 1) omega_fac(i) = omega_fac(i) + half
   if (i == 2) omega_fac(i) = omega_fac(i) - cone/(four*temp)
   if (i == 4) omega_fac(i) = omega_fac(i) + cone/(dble(48)*(temp**3))
   e_hu_migdal(:) = e_hu_migdal(:) + dble(trace_moments(:,i)*omega_fac(i))
 end do

 if (self%has_moments /= 1) then
   call destroy_matlu(matlu_tmp(:),natom)
   ABI_FREE(matlu_tmp)
 end if
 matlu_tmp => null()

 ABI_FREE(omega_fac)

 if (self%has_moments == 0) then
   call trace_prod_matlu(self_nwlo_re(:),green%occup%matlu(:),natom,trace(:))
   call destroy_matlu(self_nwlo_re(:),natom)
   ABI_FREE(self_nwlo_re)
   e_hu_migdal(:) = e_hu_migdal(:) + dble(trace(:))
 end if

 ABI_FREE(trace_moments)

 e_hu_migdal(:)  = half * e_hu_migdal(:) ! E_mig = 1/2 * Tr(Sig*G)
 e_hu_migdal_tot = sum(e_hu_migdal(:))

 ABI_FREE(trace)

 !xmig_1=zero
 !xmig_2=zero
 !xmig_3=zero

 !e_hu_migdal_tot = zero
 !do iatom=1,natom
 !  shift=czero
 !  if(paw_dmft%dmft_solv==4) shift=self%qmc_shift(iatom)+self%qmc_xmu(iatom)
!   write(std_out,*) "shiftttt",shift
  ! lpawu=paw_dmft%lpawu(iatom)
  ! if(lpawu/=-1) then
  !   xmig_1=czero
  !   xmig_2=czero
  !   xmig_3=czero
  !   ndim=2*lpawu+1
  !   do isppol=1,nsppol
  !     do ispinor = 1 , nspinor
  !       do ispinor1 = 1, nspinor
  !         do im=1,ndim
  !           do im1=1,ndim
  !             do ifreq=1,nwlo
!                write(std_out,*) ifreq,xmig_1,imag(self%oper (ifreq)%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1)),&
!&                  green%oper(ifreq)%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor )
   !              xmig_1=xmig_1 + j_dpc/beta*       &
!&                aimag(self%oper (ifreq)%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1))* &
!&                      green%oper(ifreq)%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor )* &
!&                      paw_dmft%wgt_wlo(ifreq)
!                 if(ispinor==ispinor1.and.im==im1) then
      !             se=(self%oper (ifreq)%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1)-  &
!&                      self%oper (nwlo )%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1))
!                 else
!                   se=self%oper (ifreq)%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1)
!                 endif
  !               xmig_2=xmig_2 + one/beta*real(se)* &
!&                      green%oper(ifreq)%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor )* &
!&                      paw_dmft%wgt_wlo(ifreq)
!                 if(ispinor==ispinor1.nd.im==im1.and.ifreq==1) then
  !               if(ifreq==1) then
  !                 xmig_3=xmig_3 + &
!&                   real(self%oper(nwlo )%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1)+shift)* &
!&                         green%occup%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor)/two
!                   write(std_out,*) "xmig_3",xmig_3
!                   write(std_out,*) "self",self%oper(nwlo )%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1)
!                   write(std_out,*) "shift",shift
!                   write(std_out,*) "occup", green%occup%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor)/two
  !               endif
  !             enddo
!               if(ispinor==ispinor1.and.im==im1) then
!                 xmig_3=xmig_3 + &
!&                 real(self%oper(nwlo )%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1))* &
!!&                         green%occup%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor)/two
!               endif
   !          enddo
   !        enddo
   !      enddo
   !    enddo
   !  enddo
   !  if(nsppol==1.and.nspinor==1) then
   !    e_hu_migdal(iatom)=two*real(xmig_1+xmig_2+xmig_3)
   !  else
   !    e_hu_migdal(iatom)=real(xmig_1+xmig_2+xmig_3)
   !  endif
   !  e_hu_migdal_tot = e_hu_migdal_tot + e_hu_migdal(iatom)
   !  if(abs(pawprtvol)>=3) then
   !    write(message,'(2a,3(a,5x,a,2f12.6))')ch10,&
!&         "  Interaction energy: Decomposition of Migdal energy",ch10,&
!&         "xmig_1=",xmig_1,ch10,&
!&         "xmig_2=",xmig_2,ch10,&
!&         "xmig_3=",xmig_3
!       call wrtout(std_out,message,'COLL')
 !    endif
 !  endif ! lpawu
 !enddo

end subroutine compute_migdal_energy
!!***

!!****f* m_energy/compute_dftu_energy
!! NAME
!! compute_dftu_energy
!!
!! FUNCTION
!!  Initialize noccmmp from green%occup and compute DFT+U energy with it
!!
!! INPUTS
!!  energies_dmft  = datastructure for dmft energy
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  renorm = if present change U->1 and J-> renorm just for pawuenergy
!!           renorm = J/U for the real values (does not depend on "lambda" entropy)
!!
!! OUTPUT
!!
!! SOURCE

subroutine compute_dftu_energy(energies_dmft,green,paw_dmft,pawtab,renorm)

!Arguments ------------------------------------
!type
 type(energy_type), intent(inout) :: energies_dmft
 type(green_type), intent(in) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(pawtab_type), intent(in) :: pawtab(paw_dmft%ntypat)
 real(dp), optional, intent(in) :: renorm(:)
! integer :: prtopt
!Local variables-------------------------------
 integer :: iatom,idijeff,im,im1,ims,ims1,ispinor,ispinor1,isppol,itypat
 integer :: lpawu,lpawu1,ndim,ndim1,nocc,nsploop,prt_pawuenergy
 real(dp) :: e_dc,e_dc_for_s,e_dcdc,e_dcdc_for_s,e_ee,edftumdc,edftumdc_for_s
 real(dp) :: edftumdcdc,edftumdcdc_for_s,e_ee_for_s,jpawu,upawu,xe1,xe2
 logical :: dmft_optim,t2g,x2my2d
 character(len=500) :: message
 integer, parameter :: spinor_idxs(2,4) = RESHAPE((/1,1,2,2,1,2,2,1/),(/2,4/))
 integer, parameter :: mt2g(3) = (/1,2,4/)
 real(dp), allocatable :: noccmmp(:,:,:,:),nocctot(:)
! *********************************************************************

! - allocations
! -----------------------------------------------------------------------

 e_dc       = zero
 e_dc_for_s = zero
 e_dcdc     = zero
 e_ee       = zero
 edftumdc   = zero
 edftumdcdc = zero
 nsploop    = max(paw_dmft%nsppol,paw_dmft%nspinor**2)
 nocc       = nsploop
 t2g        = (paw_dmft%dmft_t2g == 1)
 x2my2d     = (paw_dmft%dmft_x2my2d == 1)

 dmft_optim = (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7)

 isppol   = 1
 ispinor  = 1
 ispinor1 = 1

 ABI_MALLOC(nocctot,(nocc))

! - Loop and call to pawuenergy
! -----------------------------------------------------------------------
 do iatom=1,paw_dmft%natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   itypat = paw_dmft%typat(iatom)
   lpawu1 = lpawu
   if (t2g .or. x2my2d) lpawu1 = 2
   ndim  = 2*lpawu  + 1
   ndim1 = 2*lpawu1 + 1

   ABI_MALLOC(noccmmp,(2,ndim1,ndim1,nocc))
   noccmmp(:,:,:,:) = zero

! - Setup nocctot and noccmmp
! -----------------------------------------------------------------------
   nocctot(:) = zero ! contains nmmp in the n m representation
   ! Begin loop over spin/spinors to initialize noccmmp
   do idijeff=1,nsploop

     if (nsploop <= 2) then
       isppol = idijeff
     else if (nsploop == 4) then
       ispinor  = spinor_idxs(1,idijeff)
       ispinor1 = spinor_idxs(2,idijeff)
     else
       write(message,'(2a)') " BUG in m_energy: nsploop should be equal to 1, 2 or 4"
       call wrtout(std_out,message,'COLL')
     end if ! nsploop
     ! Initialize noccmmp
     do im1=1,ndim
       ims1 = im1
       ! Correct bug in computation of DFT+U energy in the t2g/x2my2d case with TRIQS
       if (x2my2d .and. dmft_optim) ims1 = 5
       if (t2g .and. dmft_optim) ims1 = mt2g(im1)
       do im=1,ndim
         ims = im
         if (x2my2d .and. dmft_optim) ims = 5
         if (t2g .and. dmft_optim) ims = mt2g(im)
         ! Here, we take the transpose in order to match pawuenergy's conventions
         noccmmp(1,ims,ims1,idijeff) = &
           & dble(green%occup%matlu(iatom)%mat(im1+(ispinor-1)*ndim,im+(ispinor1-1)*ndim,isppol))
         noccmmp(2,ims,ims1,idijeff) = &
           & aimag(green%occup%matlu(iatom)%mat(im1+(ispinor-1)*ndim,im+(ispinor1-1)*ndim,isppol))
!            noccmmp(1,im,im1,idijeff)=real(green%occup%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
!            noccmmp(2,im,im1,idijeff)=imag(green%occup%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
       end do ! im
     end do ! im1
     ! Compute nocctot
     if (green%has_charge_matlu_solver /= 2) then
       if (nsploop == 4 .and. idijeff <= 2) then
         do im1=1,ndim
           nocctot(1) = nocctot(1) + &
             & dble(green%occup%matlu(iatom)%mat(im1+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,isppol))
         end do ! im1
       else if (nsploop <= 2) then
         do im1=1,ndim
           nocctot(idijeff) = nocctot(idijeff) + &
              & dble(green%occup%matlu(iatom)%mat(im1+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,isppol))
         end do ! im1
       end if ! nsploop
     else
       if (nsploop == 4 .and. idijeff == 1) then
         nocctot(1) = green%charge_matlu_solver(2,iatom) !  total nb of elec for nspinor=2 is (iatom,2) !!
       else if (nsploop <= 2) then
         nocctot(idijeff) = green%charge_matlu_solver(idijeff,iatom) !  first spin
       end if ! nsploop
     end if ! charge_matlu_solver
   end do ! idijeff

   xe1 = e_dc
   xe2 = e_ee
    ! write(std_out,*)" nocctot(1)",nocctot(1),green%charge_matlu_solver(iatom,1)
   edftumdc   = zero
   edftumdcdc = zero
   if (present(renorm)) then
     upawu = one
     jpawu = renorm(iatom)
     prt_pawuenergy = 0
   else
     upawu = pawtab(itypat)%upawu
     jpawu = pawtab(itypat)%jpawu
     prt_pawuenergy = 3
   end if ! present(renorm)

   call pawuenergy(iatom,edftumdc,edftumdcdc,noccmmp(:,:,:,:),nocctot(:),prt_pawuenergy,pawtab(itypat),dmft_dc=paw_dmft%dmft_dc,&
                 & e_ee=e_ee,e_dc=e_dc,e_dcdc=e_dcdc,u_dmft=upawu,j_dmft=jpawu,paw_dmft=paw_dmft)

   if (paw_dmft%ientropy == 1) then
     call pawuenergy(iatom,edftumdc_for_s,edftumdcdc_for_s,noccmmp(:,:,:,:),nocctot(:),prt_pawuenergy, &
                   & pawtab(itypat),dmft_dc=paw_dmft%dmft_dc,e_ee=e_ee_for_s,e_dc=e_dc_for_s,e_dcdc=e_dcdc_for_s,&
                   & u_dmft=paw_dmft%u_for_s/Ha_eV,j_dmft=paw_dmft%j_for_s/Ha_eV)
   end if

   energies_dmft%e_dc(iatom) = e_dc - xe1
   energies_dmft%e_hu_dftu(iatom) = e_ee - xe2

   ABI_FREE(noccmmp)
 end do ! iatom

 ABI_FREE(nocctot)

! - gather results
! -----------------------------------------------------------------------
 energies_dmft%e_dc_tot = e_dc ! this is the only quantity used afterwards.
 energies_dmft%e_hu_dftu_tot = e_ee
 energies_dmft%e_dcdc = e_dcdc
 if (paw_dmft%ientropy == 1) then
   write(message,'(a,3(f14.10,3x))') "For entropy calculation E_dc_tot, u_for_s, j_for,s", &
    & e_dc_for_s,paw_dmft%u_for_s,paw_dmft%j_for_s
   call wrtout(std_out,message,'COLL')
   write(message,'(a,3(f14.10,3x))') "Reference   calculation E_dc_tot, upawu  , jpawu  ",&
      & e_dc,upawu*Ha_eV,jpawu*Ha_eV
   call wrtout(std_out,message,'COLL')
 end if ! ientropy=1

end subroutine compute_dftu_energy
!!***

!!****f* m_energy/compute_noninterentropy
!! NAME
!! compute_noninterentropy
!!
!! FUNCTION
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  green  <type(green_type)>= green function data  only for Tr(G(self-hdc))
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! SOURCE

subroutine compute_noninterentropy(cryst_struc,green,paw_dmft)

 use m_crystal, only : crystal_t

!Arguments ------------------------------------
!type
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(in) :: green
 type(paw_dmft_type), intent(inout) :: paw_dmft

!Local variables-------------------------------
 integer :: ib,ikpt,isppol,natom,nspinor,nsppol
 real(dp) :: beta,eig,fermi,s_1,s_2,occ1,occ2,f_1,e_1,f_1a,s_1a,e_2
 character(len=800) :: message
! *********************************************************************
 write(message,'(2a,i6)') ch10,"  == Compute T*Entropy for fermi level and DFT-KS eigenvalues "
 call wrtout(std_out,message,'COLL')

 natom=cryst_struc%natom
 nsppol=paw_dmft%nsppol
 nspinor=paw_dmft%nspinor
 beta=one/paw_dmft%temp
 s_1=zero
 s_1a=zero
 f_1=zero
 f_1a=zero
 e_1=zero
 e_2=zero
 s_2=zero
 do isppol=1,paw_dmft%nsppol
   do ikpt=1,paw_dmft%nkpt
     do ib=1,paw_dmft%mbandc
       eig=paw_dmft%eigen_dft(ib,ikpt,isppol)
       fermi=paw_dmft%fermie_dft
       fermi=paw_dmft%fermie
       occ1=occup_fd(eig,fermi,paw_dmft%temp)
       occ2=green%occup%ks(ib,ib,ikpt,isppol)
!       write(std_out,*) occ1,occ2

!        entropy from Fermi Dirac
       if((occ1.ge.tol9).and.((one-occ1).ge.tol9)) then
         s_1=s_1+(occ1*log(occ1)+(one-occ1)*log(one-occ1))*paw_dmft%wtk(ikpt)
!       write(std_out,*) occ1,one-occ1,"p1",(occ1*log(occ1)+(one-occ1)*log(one-occ1))*paw_dmft%wtk(ikpt)
       endif

!        Free energy from Fermi Dirac
       if((eig-fermi).ge.zero) then ! occ1 -> 0 ; 1-occ1 -> 1
         f_1=f_1-paw_dmft%wtk(ikpt)/beta*log(one+exp(-beta*(eig-fermi)))
         f_1a=f_1a+paw_dmft%wtk(ikpt)/beta*log(one-occ1)
         s_1a=s_1a+((one-occ1)*log(one-occ1)+occ1*(-beta*(eig-fermi)+log(one-occ1)))*paw_dmft%wtk(ikpt)
       else ! occ1  -> 1 , 1-occ1 -> 0
         f_1=f_1-paw_dmft%wtk(ikpt)/beta*(log(one+exp(beta*(eig-fermi)))-beta*(eig-fermi))
         f_1a=f_1a-paw_dmft%wtk(ikpt)/beta*(-log(occ1)-beta*(eig-fermi))
         s_1a=s_1a+(occ1*log(occ1)+(one-occ1)*(beta*(eig-fermi)+log(occ1)))*paw_dmft%wtk(ikpt)
       endif

!        Internal energy from Fermi Dirac
       e_1=e_1+(eig-fermi)*paw_dmft%wtk(ikpt)*occ1
       e_2=e_2+(eig-fermi)*paw_dmft%wtk(ikpt)*occ2

!        entropy from green function occupations.
       if((occ2.ge.tol9).and.((one-occ2).ge.tol9)) then
!       write(std_out,*) occ2,one-occ2,"p2",(occ2*log(occ2)+(one-occ2)*log(one-occ2))*paw_dmft%wtk(ikpt)
         s_2=s_2+(occ2*log(occ2)+(one-occ2)*log(one-occ2))*paw_dmft%wtk(ikpt)
       endif
     enddo
   enddo
 enddo
 s_1=-s_1*paw_dmft%temp
 s_1a=-s_1a*paw_dmft%temp
 s_2=-s_2*paw_dmft%temp

 write(message,'(8(2a,e20.9))') &
& ch10," T*Entropy from Fermi Dirac occupations    ", s_1,&
& ch10," T*Entropy from Fermi Dirac occupations  2 ", s_1a,&
& ch10," T*Entropy from Green function occupations ", s_2,&
& ch10," Free energy      F                        ", f_1,&
& ch10," Free energy      Fa                       ", f_1a,&
& ch10," internal energy  U                        ", e_1,&
& ch10," internal energy  U  from Gr Func Occ      ", e_2,&
& ch10," U-F                                       ", e_1-f_1
 call wrtout(std_out,message,'COLL')


end subroutine compute_noninterentropy
!!***

!!****f* m_energy/compute_free_energy
!! NAME
!! compute_free_energy
!!
!! FUNCTION
!!  Computes the different DFT+DMFT contributions to the free energy.
!!
!! INPUTS
!!  energies_dmft = datastructure for dmft energy
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  green  <type(green_type)>= green function data
!!  self  <type(self_type)>= self energy function data
!!  weiss  <type(green_type)>= weiss function data
!!  part = "band" : computes Tr(log(G_DFT)) in KS space
!!       = "main" : computes Tr(Sig*G) and Tr(log(G))
!!       = "impu" : computes the rest
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! SOURCE

subroutine compute_free_energy(energies_dmft,paw_dmft,green,part,self,weiss)

!Arguments ------------------------------------
 type(energy_type), intent(inout) :: energies_dmft
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(green_type), intent(in) :: green
 type(self_type), optional, intent(in) :: self
 type(green_type), optional, intent(in) :: weiss
 character(len=4), intent(in) :: part
!Local variables-------------------------------
 integer :: integral
 real(dp), allocatable :: e_hu_tmp(:)
! *********************************************************************

 ABI_MALLOC(e_hu_tmp,(paw_dmft%natom))

 integral = paw_dmft%dmft_triqs_compute_integral

 ! Compare Tr(log(G_DFT)) with the analytical formula (not used, simply to check that we have enough frequencies)
 if (part == "band") then
   call compute_band_energy(energies_dmft,green,paw_dmft,"nlda",fcalc_dft=1)
 end if

 if (part == "impu") then

   ! Ekin_imp
   energies_dmft%ekin_imp = green%ekin_imp

   if (integral == 1) then
     ! Tr(log(G0)) (careful, we set opt_inv to 1 since weiss contains G0^-1 rather than G0 after the dyson call)
     call compute_trace_log_loc(weiss,paw_dmft,energies_dmft%fband_weiss,opt_inv=1)
     energies_dmft%fimp = energies_dmft%fband_weiss + green%integral
   else
     energies_dmft%fimp = energies_dmft%ekin_imp + energies_dmft%e_hu_tot
   end if ! integral

   ! Integral of <dH/dlambda>
   if (integral > 0) energies_dmft%integral = green%integral

   ! Tr(log(G_imp))
   call compute_trace_log_loc(green,paw_dmft,energies_dmft%fband_imp)

   ! Tr(Sigma_imp*G_imp)
   call compute_migdal_energy(e_hu_tmp(:),energies_dmft%emig_imp,green,paw_dmft,self)
   energies_dmft%emig_imp = two * energies_dmft%emig_imp

   ! simp = entropy of the impurity
   energies_dmft%simp = (energies_dmft%ekin_imp+energies_dmft%e_hu_tot-energies_dmft%fimp) / paw_dmft%temp

 end if ! part="impu"

 if (part == "main") then

   ! Tr(Sigma*G)
   call compute_migdal_energy(e_hu_tmp(:),energies_dmft%emig_loc,green,paw_dmft,self)
   energies_dmft%emig_loc = two * energies_dmft%emig_loc

   ! Tr(log(G)) + mu*N
   energies_dmft%fband_dmft = green%trace_log

   ! fdmft = F_{dft+dmft} - E_dft
   energies_dmft%fdmft = energies_dmft%fband_dmft - energies_dmft%eband_dmft - energies_dmft%emig_loc &
      & - energies_dmft%fband_imp + energies_dmft%emig_imp + energies_dmft%fimp - energies_dmft%e_dcdc

   ! sdmft = total entropy
   energies_dmft%sdmft = (energies_dmft%edmft-energies_dmft%fdmft) / paw_dmft%temp

   call print_free_energy(energies_dmft,paw_dmft)

 end if ! part="main"

 ABI_FREE(e_hu_tmp)

end subroutine compute_free_energy
!!***

!!****f* m_energy/compute_trace_log_loc
!! NAME
!! compute_trace_log_loc
!!
!! FUNCTION
!!  Computes Tr(log(G)) in local space.
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  opt_inv = 0 (default) when green = G
!!          = 1 when green = G^-1 (useful for the Weiss field)
!!
!! OUTPUT
!!  trace = Tr(log(G_loc))
!!
!! SIDE EFFECTS
!!
!! SOURCE

subroutine compute_trace_log_loc(green,paw_dmft,trace,opt_inv)

!Arguments ------------------------------------
 type(green_type), intent(in) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 real(dp), intent(out) :: trace
 integer, optional, intent(in) :: opt_inv
!Local variables-------------------------------
 integer :: i,iatom,ierr,ifreq,info,isppol,lpawu,lwork,natom,ndim
 integer :: nmoments,nspinor,nsppol,nwlo,optinv
 real(dp) :: correction,fac,temp
 complex(dpc) :: trace_tmp
 real(dp), allocatable :: eig(:),rwork(:)
 complex(dpc), allocatable :: mat_temp(:,:),omega_fac(:),work(:)
! *********************************************************************

 optinv = 0
 if (present(opt_inv)) optinv = opt_inv

 natom    = paw_dmft%natom
 nmoments = green%nmoments - 1
 nspinor  = paw_dmft%nspinor
 nsppol   = paw_dmft%nsppol
 nwlo     = green%nw
 temp     = paw_dmft%temp
 trace    = zero
 ndim     = nspinor * (2*paw_dmft%maxlpawu+1)

 ABI_MALLOC(eig,(ndim))
 ABI_MALLOC(rwork,(3*ndim-2))
 ABI_MALLOC(work,(2*ndim-1))
 ABI_MALLOC(mat_temp,(ndim,ndim))
 call zheev('n','u',ndim,mat_temp(:,:),ndim,eig(:),work(:),-1,rwork(:),info)
 lwork = int(work(1))
 ABI_FREE(work)
 ABI_MALLOC(work,(lwork))
 ABI_FREE(mat_temp)

 do ifreq=1,nwlo
   if (green%distrib%procf(ifreq) /= paw_dmft%myproc) cycle
   fac = merge(temp*two,temp,nsppol==1.and.nspinor==1)
   trace_tmp = czero
   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     ndim = nspinor * (2*lpawu+1)
     ABI_MALLOC(mat_temp,(ndim,ndim))
     do isppol=1,nsppol
       call abi_xgemm("n","c",ndim,ndim,ndim,cone,green%oper(ifreq)%matlu(iatom)%mat(:,:,isppol),ndim,&
                    & green%oper(ifreq)%matlu(iatom)%mat(:,:,isppol),ndim,czero,mat_temp(:,:),ndim)
       call zheev('n','u',ndim,mat_temp(:,:),ndim,eig(:),work(:),lwork,rwork(1:3*ndim-2),info)

       if (optinv == 1) then
         trace_tmp = trace_tmp - sum(log(eig(1:ndim)/(paw_dmft%omega_lo(ifreq)**2)))
       else
         trace_tmp = trace_tmp + sum(log(eig(1:ndim)*(paw_dmft%omega_lo(ifreq)**2)))
       end if

     end do ! isppol
     if (ifreq == nwlo) then
       correction = fac * nsppol * ndim * log(two)
       trace = trace - correction
     end if
     ABI_FREE(mat_temp)
   end do ! iatom
   trace = trace + dble(trace_tmp)*fac
 end do ! ifreq

 ABI_FREE(rwork)
 ABI_FREE(work)
 ABI_FREE(eig)

 call xmpi_sum(trace,paw_dmft%spacecomm,ierr)

 ABI_MALLOC(omega_fac,(nmoments))

 do i=1,nmoments
   omega_fac(i) = czero
   do ifreq=nwlo,1,-1 ! NEVER change the summation order and DON'T use the intrinsic SUM
     omega_fac(i) = omega_fac(i) + cone / (paw_dmft%omega_lo(ifreq))**i
   end do
   omega_fac(i) = - two * temp * omega_fac(i) / (j_dpc)**i
   if (i == 1) omega_fac(i) = omega_fac(i) + half
   if (i == 2) omega_fac(i) = omega_fac(i) - cone/(four*temp)
   if (i == 4) omega_fac(i) = omega_fac(i) + cone/(dble(48)*(temp**3))
 end do ! i

 ! Do not use dot_product
 trace = trace + dble(sum(green%trace_moments_log_loc(1:nmoments)*omega_fac(1:nmoments)))

 ABI_FREE(omega_fac)

end subroutine compute_trace_log_loc
!!***

!!****f* m_energy/print_free_energy
!! NAME
!! print_free_energy
!!
!! FUNCTION
!!  Prints the different DFT+DMFT contributions to the free energy.
!!
!! INPUTS
!!  energies_dmft = datastructure for dmft energy
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! SOURCE

subroutine print_free_energy(energies_dmft,paw_dmft)

!Arguments ------------------------------------
 type(energy_type), intent(in) :: energies_dmft
 type(paw_dmft_type), intent(in) :: paw_dmft
!Local variables-------------------------------
 integer :: integral
 real(dp) :: temp
 character(len=10000) :: message,message2
! *********************************************************************

 integral = paw_dmft%dmft_triqs_compute_integral
 temp = paw_dmft%temp

 write(message,'(a,5x,2a,5x,a,10(a,5x,a,2x,f18.11),a)') ch10, &
     & "-----------------------------------------------",ch10, &
     & "--- Free Energy in DMFT (in Ha)  ",ch10, &
     & "--- E_hu                      (1) (Ha.) = ",energies_dmft%e_hu_tot,ch10, &
     & "--- E_dc                      (2) (Ha.) = ",energies_dmft%e_dc_tot,ch10, &
     & "--- E_dmft                (1)-(2) (Ha.) = ",energies_dmft%edmft,ch10, &
     & "--- Tr(log(G))+mu*N           (3) (Ha.) = ",energies_dmft%fband_dmft,ch10, &
     & "--- Tr(Sigma*G)               (4) (Ha.) = ",energies_dmft%emig_loc,ch10, &
     & "--- Tr(V_dc*rho)              (5) (Ha.) = ",-energies_dmft%e_dcdc+energies_dmft%e_dc_tot,ch10, &
     & "--- E_band_dmft               (6) (Ha.) = ",energies_dmft%eband_dmft,ch10, &
     & "--- Tr(log(Gimp))             (7) (Ha.) = ",energies_dmft%fband_imp,ch10, &
     & "--- Tr(Sigma_imp*G_imp)       (8) (Ha.) = ",energies_dmft%emig_imp,ch10, &
     & "--- E_kinetic_imp             (9) (Ha.) = ",energies_dmft%ekin_imp,ch10

 if (integral > 0) then
   write(message2,'(8(5x,a,2x,f18.11,a),5x,a)') &
     & "--- Tr(log(G0))              (10) (Ha.) = ",energies_dmft%fband_weiss,ch10, &
     & "--- Integral                 (11) (Ha.) = ",energies_dmft%integral,ch10, &
     & "--- F_imp               (10)+(11) (Ha.) = ",energies_dmft%fimp,ch10, &
     & "--- (-kT)*S_imp (10)+(11)-(9)-(1) (Ha.) = ",-temp*energies_dmft%simp,ch10, &
     & "--- S_imp                    (12)       = ",energies_dmft%simp,ch10, &
     & "--- F_dmft                   (13) (Ha.) = ",energies_dmft%fdmft,ch10, &
     & "--- (-kT)*S_dmft     (13)-(1)+(2) (Ha.) = ",-temp*energies_dmft%sdmft,ch10, &
     & "--- S_dmft                              = ",energies_dmft%sdmft,ch10, &
     & "-----------------------------------------------"
 else
   write(message2,'(3(5x,a,2x,f18.11,a),5x,a)') &
     & "--- F_dmft+T*S_imp           (10) (Ha.) = ",energies_dmft%fdmft,ch10, &
     & "--- (-kT)*(S_dmft-S_imp)     (11) (Ha.) = ",-temp*energies_dmft%sdmft,ch10, &
     & "--- S_dmft-S_imp             (12)       = ",energies_dmft%sdmft,ch10, &
     & "-----------------------------------------------"
 end if ! integral

 call wrtout(std_out,trim(adjustl(message))//trim(message2),'COLL')

end subroutine print_free_energy
!!***

END MODULE m_energy
!!***

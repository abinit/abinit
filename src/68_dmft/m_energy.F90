!!****m* ABINIT/m_energy
!! NAME
!!  m_energy
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2020 ABINIT group (BAmadon)
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

MODULE m_energy

 use defs_basis
 use m_errors
 use m_abicore

 use m_pawtab, only : pawtab_type
 use m_paw_correlations, only : pawuenergy

 use m_green, only : green_type,icip_green,destroy_green,compa_occup_ks
 use m_self, only : self_type,initialize_self,destroy_self,print_self,new_self,make_qmcshift_self
 use m_paw_dmft, only : paw_dmft_type

 use m_matlu, only : matlu_type,init_matlu,prod_matlu,diag_matlu,destroy_matlu,conjg_matlu,&
& ln_matlu,add_matlu,zero_matlu,shift_matlu,copy_matlu,trace_matlu,print_matlu
 use m_oper, only : trace_oper
 use m_crystal, only : crystal_t

 implicit none

 private

 public :: init_energy
 public :: compute_energy
 public :: compute_migdal_energy
 public :: compute_dftu_energy
 public :: destroy_energy
 public :: print_energy
 public :: compute_noninterentropy

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

  real(dp) :: eband_dft

  real(dp) :: eband_dmft

  real(dp) :: e_dc_tot

  real(dp) :: e_hu_tot

  real(dp) :: e_hu_dftu_tot

  real(dp) :: e_hu_mig_tot

  real(dp) :: e_hu_qmc_tot

  real(dp) :: edmft

  real(dp) :: natom

  real(dp), allocatable :: e_dc(:)

  real(dp), allocatable :: e_hu(:)

  real(dp), allocatable :: e_hu_dftu(:)

  real(dp), allocatable :: e_hu_mig(:)

  real(dp), allocatable :: e_hu_qmc(:)

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
!! energies_dmft  = structure of data for dmft of type energy_type
!!
!! PARENTS
!!      m_dmft
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine init_energy(cryst_struc,energies_dmft)

!Arguments ------------------------------------
!type
 type(crystal_t),intent(in) :: cryst_struc
 type(energy_type), intent(inout) :: energies_dmft
!Local variables ------------------------------------
!************************************************************************

 ABI_ALLOCATE(energies_dmft%e_dc,(cryst_struc%natom))
 ABI_ALLOCATE(energies_dmft%e_hu,(cryst_struc%natom))
 ABI_ALLOCATE(energies_dmft%e_hu_dftu,(cryst_struc%natom))
 ABI_ALLOCATE(energies_dmft%e_hu_mig,(cryst_struc%natom))
 ABI_ALLOCATE(energies_dmft%e_hu_qmc,(cryst_struc%natom))
 energies_dmft%e_dc=zero
 energies_dmft%e_hu=zero
 energies_dmft%e_hu_dftu=zero
 energies_dmft%e_hu_mig=zero
 energies_dmft%e_hu_qmc=zero
 energies_dmft%eband_dft=zero
 energies_dmft%eband_dmft=zero
 energies_dmft%e_dc_tot=zero
 energies_dmft%e_hu_tot=zero
 energies_dmft%e_hu_dftu_tot=zero
 energies_dmft%e_hu_mig_tot=zero
 energies_dmft%e_hu_qmc_tot=zero
 energies_dmft%edmft=zero
 energies_dmft%natom=cryst_struc%natom

end subroutine init_energy
!!***

!!****f* m_energy/destroy_energy
!! NAME
!! destroy_energy
!!
!! FUNCTION
!!  deallocate energies_dmft
!!
!! INPUTS
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dmft
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine destroy_energy(energies_dmft,paw_dmft)

!Arguments ------------------------------------
!scalars
 type(energy_type),intent(inout) :: energies_dmft
 type(paw_dmft_type), intent(inout) :: paw_dmft
!Local variables-------------------------------
! *********************************************************************
  paw_dmft%edmft=energies_dmft%edmft
 if ( allocated(energies_dmft%e_dc) )   then
   ABI_DEALLOCATE(energies_dmft%e_dc)
 end if
 if ( allocated(energies_dmft%e_hu) )   then
   ABI_DEALLOCATE(energies_dmft%e_hu)
 end if
 if ( allocated(energies_dmft%e_hu_dftu) )  then
   ABI_DEALLOCATE(energies_dmft%e_hu_dftu)
 end if
 if ( allocated(energies_dmft%e_hu_mig) )  then
   ABI_DEALLOCATE(energies_dmft%e_hu_mig)
 end if
  if ( allocated(energies_dmft%e_hu_qmc) )  then
   ABI_DEALLOCATE(energies_dmft%e_hu_qmc)
 end if


end subroutine destroy_energy
!!***

!!****f* m_energy/print_energy
!! NAME
!! print_energy
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_energy
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine print_energy(cryst_struc,energies_dmft,pawprtvol,pawtab,idmftloop)

!Arguments ------------------------------------
!type
 type(crystal_t),intent(in) :: cryst_struc
 type(energy_type),intent(in) :: energies_dmft
 type(pawtab_type),intent(in)  :: pawtab(cryst_struc%ntypat)
 integer, intent(in) :: pawprtvol,idmftloop

!Local variables-------------------------------
 integer :: iatom
 character(len=1000) :: message
! *********************************************************************
 if(abs(pawprtvol)>=3) then
   do iatom=1,cryst_struc%natom
     if(pawtab(cryst_struc%typat(iatom))%lpawu/=-1) then
       write(message,'(a,4x,a,i3,a,5x,f12.6)')  &
&       ch10,"For Correlated Atom",iatom,", E_hu =",energies_dmft%e_hu(iatom)
       call wrtout(std_out,message,'COLL')
       write(message,'(26x,a,1x,f12.6)')  &
&       ", E_hu_mig =",energies_dmft%e_hu_mig(iatom)
       call wrtout(std_out,message,'COLL')
       write(message,'(26x,a,f12.6)')  &
&       ", E_hu_qmc  =",energies_dmft%e_hu_qmc(iatom)
       call wrtout(std_out,message,'COLL')
       write(message,'(26x,a,f12.6)')  &
&       ", E_hu_dftu =",energies_dmft%e_hu_dftu(iatom)
       call wrtout(std_out,message,'COLL')
       write(message,'(26x,a,f12.6)')  &
&       ", E_dc =",energies_dmft%e_dc(iatom)
       call wrtout(std_out,message,'COLL')
     endif
   enddo
 endif
 write(message,'(a,5x,2a,5x,a,9(a,5x,a,2x,f15.11),a,5x,a)') ch10 &
&      ,"-----------------------------------------------",ch10 &
&      ,"--- Energy in DMFT (in Ha)  ",ch10 &
&      ,"--- E_bandlda (1)  (Ha.) = ",energies_dmft%eband_dft,ch10 &
&      ,"--- E_banddmft(2)  (Ha.) = ",energies_dmft%eband_dmft,ch10 &
&      ,"--- E_hu      (3)  (Ha.) = ",energies_dmft%e_hu_tot,ch10 &
&      ,"--- E_hu_mig  (4)  (Ha.) = ",energies_dmft%e_hu_mig_tot,ch10 &
&      ,"--- E_hu_qmc  (4)  (Ha.) = ",energies_dmft%e_hu_qmc_tot,ch10 &
&      ,"--- E_hu_dftu (5)  (Ha.) = ",energies_dmft%e_hu_dftu_tot,ch10 &
&      ,"--- E_dc      (6)  (Ha.) = ",energies_dmft%e_dc_tot,ch10 &
&      ,"--- edmft=(    3-6)(Ha.) = ",energies_dmft%edmft,ch10 &
&      ,"---       (2-1+3-6)(Ha.) = ",energies_dmft%eband_dmft-energies_dmft%eband_dft+energies_dmft%edmft,ch10 &
&      ,"-----------------------------------------------"
 call wrtout(std_out,message,'COLL')
 if(idmftloop>=1) then
   write(message,'(a,i3,1x,f15.11,a)') " (Edmft",idmftloop,energies_dmft%edmft,")"
   call wrtout(ab_out,message,'COLL')
 endif
end subroutine print_energy
!!***

!!****f* m_energy/compute_energy
!! NAME
!! compute_energy
!!
!! FUNCTION
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  energies_dmft <type(energy_type)> = DMFT energy structure data
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawprtvol
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  self  <type(self_type)>= self energy function data
!!  occ_type=  character ("lda" or "nlda") for printing.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! energies_dmft <type(energy_type)> = DMFT energy structure data
!!
!! PARENTS
!!      m_dmft
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine compute_energy(cryst_struc,energies_dmft,green,paw_dmft,pawprtvol,pawtab,self,occ_type,part)

!Arguments ------------------------------------
!type
 type(energy_type),intent(inout) :: energies_dmft
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(in) :: green
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(pawtab_type),intent(in)  :: pawtab(cryst_struc%ntypat)
 type(self_type), intent(inout) :: self
 integer, intent(in) :: pawprtvol
 character(len=4), intent(in) :: occ_type
 character(len=4), intent(in) :: part
! integer :: prtopt

!Local variables-------------------------------
 integer :: iatom,lpawu
 integer :: natom,nspinor,nsppol
 real(dp) :: beta
 character(len=500) :: message
 real(dp), allocatable :: e_hu_migdal(:)
 real(dp) :: e_hu_migdal_tot
! *********************************************************************
 if(part=='both') then
   write(message,'(2a)') ch10,"  == Compute DFT+DMFT energy terms "
   call wrtout(std_out,message,'COLL')
 else if(part=='band') then
   write(message,'(2a)') ch10,"  == Compute DFT+DMFT energy terms : Band energy terms"
   call wrtout(std_out,message,'COLL')
 else if(part=='corr') then
   write(message,'(2a)') ch10,"  == Compute DFT+DMFT energy terms : Correlation energy terms only"
   call wrtout(std_out,message,'COLL')
 else if(part=='none') then
 endif

! Only imaginary frequencies here
 if(green%w_type=="real".or.self%w_type=="real") then
   message = 'compute_energy not implemented for real frequency'
   MSG_BUG(message)
 endif
 natom=cryst_struc%natom
 nsppol=paw_dmft%nsppol
 nspinor=paw_dmft%nspinor
 beta=one/paw_dmft%temp

 if(part=='band'.or.part=='both') then
! == Compute Band Energy Alternative version: two steps
!                 == Compute Tr[ln G^{-1}] and -Tr[(Self-hdc)G_dmft]
! -----------------------------------------------------------------------
   if(part=='band') then ! ie if thdyn="fcalc" in m_dmft.F90
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
     call compute_band_energy(energies_dmft,green,paw_dmft,occ_type,ecalc_dft=1)
!     write(message,'(2a,f10.6)') ch10,"Compute Band energy ref  -ef          ",energies_dmft%eband_dmft
!     call wrtout(std_out,message,'COLL')
!! call wrtout(ab_out,message,'COLL')
!     write(message,'(2a,f10.6)') ch10,"Compute Band energy ref   lda         ",energies_dmft%eband_dft
!     call wrtout(std_out,message,'COLL')
!! if(occ_type=="nlda") eband2=energies_dmft%eband_dmft
   else
     call compute_band_energy(energies_dmft,green,paw_dmft,occ_type,ecalc_dft=1)
   endif

 endif
 if(part=='corr'.or.part=='both') then

! == Compute Correlation energy from Migdal formula
! -----------------------------------------------------------------------
   ABI_ALLOCATE(e_hu_migdal,(cryst_struc%natom))
   e_hu_migdal(:) = zero
   call compute_migdal_energy(cryst_struc,e_hu_migdal,e_hu_migdal_tot,green,paw_dmft,pawprtvol,self)
   energies_dmft%e_hu_mig(:)= e_hu_migdal(:)
   energies_dmft%e_hu_mig_tot = e_hu_migdal_tot
! write(std_out,*) "MIGDAL",e_hu_migdal_tot,e_hu_migdal
   ABI_DEALLOCATE(e_hu_migdal)

! == Compute Correlation energy from QMC correlations.
! -----------------------------------------------------------------------
   energies_dmft%e_hu_qmc_tot = zero
   do iatom=1,natom
     lpawu=paw_dmft%lpawu(iatom)
     if(lpawu/=-1) then
       if(paw_dmft%dmft_solv==4.or.paw_dmft%dmft_solv==5.or.paw_dmft%dmft_solv==8) then
         energies_dmft%e_hu_qmc(iatom) = green%ecorr_qmc(iatom)
         energies_dmft%e_hu_qmc_tot = energies_dmft%e_hu_qmc_tot + green%ecorr_qmc(iatom)
       endif
     endif ! lpawu
   enddo ! iatom

! == Compute DFT+U interaction energy
! -----------------------------------------------------------------------
   call compute_dftu_energy(cryst_struc,energies_dmft,green,paw_dmft,pawtab)
   if(abs(paw_dmft%dmft_solv)<=1) then
     energies_dmft%e_hu= energies_dmft%e_hu_dftu
     energies_dmft%e_hu_tot= energies_dmft%e_hu_dftu_tot
     if((abs(energies_dmft%e_hu_tot-energies_dmft%e_hu_mig_tot).ge.0.000001).and.(occ_type/=" lda")) then
       write(message,'(2a,2e18.8,a)') ch10,'   BUG: Migdal energy and DFT+U energy do not coincide',&
&       energies_dmft%e_hu_tot,energies_dmft%e_hu_mig_tot,occ_type
       MSG_ERROR(message)
     endif
   else if(paw_dmft%dmft_solv==2.or.paw_dmft%dmft_solv==6.or.paw_dmft%dmft_solv==7.or.paw_dmft%dmft_solv==9) then
     energies_dmft%e_hu= energies_dmft%e_hu_mig
     energies_dmft%e_hu_tot= energies_dmft%e_hu_mig_tot
     energies_dmft%e_hu_qmc_tot = energies_dmft%e_hu_tot
   else if(paw_dmft%dmft_solv==4.or.paw_dmft%dmft_solv==5.or.paw_dmft%dmft_solv==8) then
     if(paw_dmft%dmft_solv==8) then
       write(message,'(2a)') ch10,"Warning, energy is recently computed, not checked"
       call wrtout(std_out,message,'COLL')
     endif
     energies_dmft%e_hu= energies_dmft%e_hu_qmc
     energies_dmft%e_hu_tot= energies_dmft%e_hu_qmc_tot
   endif
!   energies_dmft%edmft=energies_dmft%e_hu_mig_tot-energies_dmft%e_dc_tot
   energies_dmft%edmft=energies_dmft%e_hu_tot-energies_dmft%e_dc_tot


 endif ! part

! if(part='corr'.or.part='both') then
 if(part/='none') then
   call print_energy(cryst_struc,energies_dmft,pawprtvol,pawtab,paw_dmft%idmftloop)
 endif
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
!!
!! INPUTS
!!  energies_dmft <type(energy_type)> = DMFT energy structure data
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  occ_type=  character ("lda" or "nlda") for printing.
!!  fcalc_dft= if present, compute free energy instead of total energy.
!!  fcalc_dft= 2 do a free energy calculation with the dmft fermie level
!!           = 2/3 subtract fermie to eigenvalues and free energy
!!           =1/3 fermie_dft level and free energy
!!  ecalc_dft= 1 subtract fermie to eigenvalues
!! OUTPUT
!!
!! SIDE EFFECTS
!!  energies_dmft <type(energy_type)> = DMFT energy structure data
!!
!! PARENTS
!!      m_energy
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine compute_band_energy(energies_dmft,green,paw_dmft,occ_type,ecalc_dft,fcalc_dft,ecalc_dmft)

!Arguments ------------------------------------
!type
 type(energy_type),intent(inout) :: energies_dmft
 type(green_type),intent(in) :: green
 type(paw_dmft_type), intent(inout) :: paw_dmft
 character(len=4), intent(in) :: occ_type
 integer, intent(in), optional :: ecalc_dft
 integer, intent(in), optional :: fcalc_dft
 integer, intent(in), optional :: ecalc_dmft
! integer :: prtopt

!Local variables-------------------------------
 integer :: ib,ikpt,isppol
 real(dp) :: beta , fermie_used, totch,totch2,totch3
 character(len=500) :: message
! *********************************************************************
 write(message,'(2a)') ch10,"  == Compute Band Energy terms for DMFT "
 call wrtout(std_out,message,'COLL')
 beta=one/paw_dmft%temp

! == Compute Band Energy
! -----------------------------------------------------------------------
 energies_dmft%eband_dft=zero
 energies_dmft%eband_dmft=zero
 totch=zero
 totch2=zero
 totch3=zero
 do isppol=1,paw_dmft%nsppol
   do ikpt=1,paw_dmft%nkpt
     do ib=1,paw_dmft%mbandc
       if(present(fcalc_dft)) then
         if (fcalc_dft==1.or.fcalc_dft==3) fermie_used=paw_dmft%fermie_dft
         if (fcalc_dft==2.or.fcalc_dft==4) fermie_used=paw_dmft%fermie ! only for B3 terms
         if((paw_dmft%eigen_dft(isppol,ikpt,ib)-fermie_used).ge.zero) then
           energies_dmft%eband_dft=energies_dmft%eband_dft - &
&             log ( one + exp ( - beta* (paw_dmft%eigen_dft(isppol,ikpt,ib)-fermie_used)))*paw_dmft%wtk(ikpt)
         else
           energies_dmft%eband_dft=energies_dmft%eband_dft  &
&             -(log ( one + exp ( beta* (paw_dmft%eigen_dft(isppol,ikpt,ib)-fermie_used))))*paw_dmft%wtk(ikpt) &
&              + beta*paw_dmft%eigen_dft(isppol,ikpt,ib)*paw_dmft%wtk(ikpt)
           if(fcalc_dft==3.or.fcalc_dft==2) then
                   ! subtract fermi level, (useful to directly count the number of electrons)
             energies_dmft%eband_dft=energies_dmft%eband_dft  &
&                                 - beta*fermie_used*paw_dmft%wtk(ikpt)
             totch=totch+paw_dmft%wtk(ikpt)
           endif
         endif
       else ! usual calculation: total non interacting energy
         fermie_used=paw_dmft%fermie_dft
!            write(std_out,*) "isppol,ikpt,ib",isppol,ikpt,ib
!            write(std_out,*) "paw_dmft%eigen_dft",paw_dmft%eigen_dft(isppol,ikpt,ib)
!            write(std_out,*) green%occup%ks(isppol,ikpt,ib,ib)
!            write(std_out,*) occup_fd(paw_dmft%eigen_dft(isppol,ikpt,ib),paw_dmft%fermie,paw_dmft%temp)
         if(present(ecalc_dft)) then
           if(ecalc_dft==1.or.ecalc_dft==3) fermie_used=paw_dmft%fermie_dft
           if(ecalc_dft==2.or.ecalc_dft==4) fermie_used=paw_dmft%fermie ! only for B3 terms
           if(ecalc_dft==3.or.ecalc_dft==2) then
             energies_dmft%eband_dft=energies_dmft%eband_dft- &
&               occup_fd(paw_dmft%eigen_dft(isppol,ikpt,ib),fermie_used,paw_dmft%temp)*&
&               fermie_used*paw_dmft%wtk(ikpt)
             totch2=totch2+paw_dmft%wtk(ikpt)*occup_fd(paw_dmft%eigen_dft(isppol,ikpt,ib),fermie_used,paw_dmft%temp)
           endif
         endif
         energies_dmft%eband_dft=energies_dmft%eband_dft+ &
&           occup_fd(paw_dmft%eigen_dft(isppol,ikpt,ib),fermie_used,paw_dmft%temp)*&
&           paw_dmft%eigen_dft(isppol,ikpt,ib)*paw_dmft%wtk(ikpt)
       endif
       energies_dmft%eband_dmft=energies_dmft%eband_dmft+ &
&         green%occup%ks(isppol,ikpt,ib,ib)*&
&         paw_dmft%eigen_dft(isppol,ikpt,ib)*paw_dmft%wtk(ikpt)
         totch3=totch3+paw_dmft%wtk(ikpt)*green%occup%ks(isppol,ikpt,ib,ib)
         if(present(ecalc_dmft)) then
           energies_dmft%eband_dmft=energies_dmft%eband_dmft- &
&              green%occup%ks(isppol,ikpt,ib,ib)*&
&              paw_dmft%fermie*paw_dmft%wtk(ikpt)
         endif
     enddo
   enddo
 enddo
 if(paw_dmft%nsppol==1.and.paw_dmft%nspinor==1)  energies_dmft%eband_dft=two*energies_dmft%eband_dft
 if(paw_dmft%nsppol==1.and.paw_dmft%nspinor==1)  energies_dmft%eband_dmft=two*energies_dmft%eband_dmft
 if(present(fcalc_dft)) then
   energies_dmft%eband_dft=energies_dmft%eband_dft/beta
   if(fcalc_dft==3.or.fcalc_dft==2) write(std_out,*) "compute_band_energy totch",totch
 endif
 if(present(ecalc_dft)) then
   if(ecalc_dft==3.or.ecalc_dft==2) write(std_out,*) "compute_band_energy totch2",totch2
 endif
! write(std_out,*) "compute_band_energy totch3",totch3

 if (occ_type==" lda") then
   if(abs(energies_dmft%eband_dft-energies_dmft%eband_dmft)>tol5) then
     write(message,'(5x,a,a,a,15x,a,f12.6,a,15x,a,5x,f12.5)')  "Warning !:"&
&     ,"Differences between band energy from DFT occupations",ch10&
&     ,"and DFT green function is:",energies_dmft%eband_dft-energies_dmft%eband_dmft,ch10&
&     ,"which is larger than",tol5
     call wrtout(std_out,message,'COLL')
     write(message,'(a)') &
&     "   Action: increase number of frequencies, or reduce the number of high energies_dmft bands"
     call wrtout(std_out,message,'COLL')
   else
     write(message,'(a,a,a,10x,a,f12.6,a,10x,a,5x,f12.5)')  "          "&
&     ,"Differences between band energy from DFT occupations",ch10&
&     ,"and DFT green function is:",energies_dmft%eband_dft-energies_dmft%eband_dmft,ch10&
&     ,"which is smaller than",tol5
     call wrtout(std_out,message,'COLL')
   endif
 endif


end subroutine compute_band_energy
!!***

!!****f* m_energy/compute_migdal_energy
!! NAME
!! compute_migdal_energy
!!
!! FUNCTION
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)> = crystal structure data
!!  energies_dmft <type(energy_type)> = DMFT energy structure data
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawprtvol
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  self  <type(self_type)>= self energy function data
!!
!! OUTPUT
!!  e_hu_mig(natom)= Migdal energy for each atom.
!!  e_hu_mig_tot= Total Migdal energy.
!!
!! PARENTS
!!      m_energy
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine compute_migdal_energy(cryst_struc,e_hu_migdal,e_hu_migdal_tot,green,paw_dmft,pawprtvol,self)

#ifdef FC_INTEL
!DEC$ NOOPTIMIZE
#endif

!Arguments ------------------------------------
!type
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(in) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 real(dp),intent(out) :: e_hu_migdal_tot
 real(dp),intent(out) :: e_hu_migdal(paw_dmft%natom)
 type(self_type), intent(in) :: self
 integer, intent(in) :: pawprtvol
! integer :: prtopt

!Local variables-------------------------------
 integer :: iatom,ifreq,im,im1,ispinor,ispinor1,isppol,lpawu
 integer :: natom,ndim,nspinor,nsppol,nwlo
 real(dp) :: beta
 complex(dpc) :: xmig_1,xmig_2,xmig_3,se,shift
 character(len=500) :: message
! *********************************************************************

! Only imaginary frequencies here
 if(green%w_type=="real".or.self%w_type=="real") then
   message = 'compute_migdal_energy not implemented for real frequency'
   MSG_BUG(message)
 endif

! == Compute Correlation energy from Migdal formula
! -----------------------------------------------------------------------
 natom=cryst_struc%natom
 nsppol=paw_dmft%nsppol
 nspinor=paw_dmft%nspinor
 beta=one/paw_dmft%temp
 nwlo=green%nw
 if (green%nw/=self%nw) then
   message = 'self and green do not contains the same number of frequencies'
   MSG_BUG(message)
 endif
! write(std_out,*) "beta",beta

 xmig_1=zero
 xmig_2=zero
 xmig_3=zero

 e_hu_migdal_tot = zero
 do iatom=1,natom
   shift=czero
   if(paw_dmft%dmft_solv==4) shift=self%qmc_shift(iatom)+self%qmc_xmu(iatom)
!   write(std_out,*) "shiftttt",shift
   lpawu=paw_dmft%lpawu(iatom)
   if(lpawu/=-1) then
     xmig_1=czero
     xmig_2=czero
     xmig_3=czero
     ndim=2*lpawu+1
     do isppol=1,nsppol
       do ispinor = 1 , nspinor
         do ispinor1 = 1, nspinor
           do im=1,ndim
             do im1=1,ndim
               do ifreq=1,nwlo
!                write(std_out,*) ifreq,xmig_1,imag(self%oper (ifreq)%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1)),&
!&                  green%oper(ifreq)%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor )
                 xmig_1=xmig_1 + j_dpc/beta*       &
&                aimag(self%oper (ifreq)%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1))* &
&                      green%oper(ifreq)%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor )* &
&                      paw_dmft%wgt_wlo(ifreq)
!                 if(ispinor==ispinor1.and.im==im1) then
                   se=(self%oper (ifreq)%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1)-  &
&                      self%oper (nwlo )%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1))
!                 else
!                   se=self%oper (ifreq)%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1)
!                 endif
                 xmig_2=xmig_2 + one/beta*real(se)* &
&                      green%oper(ifreq)%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor )* &
&                      paw_dmft%wgt_wlo(ifreq)
!                 if(ispinor==ispinor1.nd.im==im1.and.ifreq==1) then
                 if(ifreq==1) then
                   xmig_3=xmig_3 + &
&                   real(self%oper(nwlo )%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1)+shift)* &
&                         green%occup%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor)/two
!                   write(std_out,*) "xmig_3",xmig_3
!                   write(std_out,*) "self",self%oper(nwlo )%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1)
!                   write(std_out,*) "shift",shift
!                   write(std_out,*) "occup", green%occup%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor)/two
                 endif
               enddo
!               if(ispinor==ispinor1.and.im==im1) then
!                 xmig_3=xmig_3 + &
!&                 real(self%oper(nwlo )%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1))* &
!!&                         green%occup%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor)/two
!               endif
             enddo
           enddo
         enddo
       enddo
     enddo
     if(nsppol==1.and.nspinor==1) then
       e_hu_migdal(iatom)=two*real(xmig_1+xmig_2+xmig_3)
     else
       e_hu_migdal(iatom)=real(xmig_1+xmig_2+xmig_3)
     endif
     e_hu_migdal_tot = e_hu_migdal_tot + e_hu_migdal(iatom)
     if(abs(pawprtvol)>=3) then
       write(message,'(2a,3(a,5x,a,2f12.6))')ch10,&
&         "  Interaction energy: Decomposition of Migdal energy",ch10,&
&         "xmig_1=",xmig_1,ch10,&
&         "xmig_2=",xmig_2,ch10,&
&         "xmig_3=",xmig_3
       call wrtout(std_out,message,'COLL')
     endif
   endif ! lpawu
 enddo

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
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  renorm = if present change U->1 and J-> renorm just for pawuenergy
!!           renorm = J/U for the real values (does not depend on "lambda" entropy)
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dmft,m_energy
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine compute_dftu_energy(cryst_struc,energies_dmft,green,paw_dmft,pawtab,renorm)

!Arguments ------------------------------------
!type
 type(energy_type),intent(inout) :: energies_dmft
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(in) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(pawtab_type),target,intent(in)  :: pawtab(cryst_struc%ntypat)
 real(dp), optional, intent(in) :: renorm(:)
! integer :: prtopt

!Local variables-------------------------------
 integer :: iatom,idijeff,im,im1,ispinor,ispinor1,isppol,ldim,lpawu
 integer :: nocc,nsploop,prt_pawuenergy
 real(dp) :: upawu,jpawu
 real(dp) :: edftumdcdc,edftumdc,e_ee,e_dc,e_dcdc,xe1,xe2
 character(len=500) :: message
! arrays
 integer,parameter :: spinor_idxs(2,4)=RESHAPE((/1,1,2,2,1,2,2,1/),(/2,4/))
 real(dp),allocatable :: noccmmp(:,:,:,:),nocctot(:)
 type(pawtab_type),pointer :: pawtab_

! *********************************************************************

! - allocations
! -----------------------------------------------------------------------

 nsploop=max(paw_dmft%nsppol,paw_dmft%nspinor**2)
 nocc=nsploop
 e_ee=zero
 e_dc=zero
 e_dcdc=zero
 edftumdc=zero
 edftumdcdc=zero
 isppol=0
 ispinor=0
 ispinor1=0

! - Loop and call to pawuenergy
! -----------------------------------------------------------------------
 do iatom=1,cryst_struc%natom
   pawtab_ => pawtab(cryst_struc%typat(iatom))
   lpawu=paw_dmft%lpawu(iatom)
   if(lpawu.ne.-1) then
     ldim=2*lpawu+1

     ABI_ALLOCATE(noccmmp,(2,2*pawtab_%lpawu+1,2*pawtab_%lpawu+1,nocc))
     ABI_ALLOCATE(nocctot,(nocc))
     noccmmp(:,:,:,:)=zero ; nocctot(:)=zero

! - Setup nocctot and noccmmp
! -----------------------------------------------------------------------
     nocctot(:)=zero ! contains nmmp in the n m representation
! Begin loop over spin/spinors to initialize noccmmp
     do idijeff=1,nsploop
       if(nsploop==2) then
         isppol=spinor_idxs(1,idijeff)
         ispinor=1
         ispinor1=1
       else if(nsploop==4) then
         isppol=1
         ispinor=spinor_idxs(1,idijeff)
         ispinor1=spinor_idxs(2,idijeff)
       else if(nsploop==1) then
         isppol=1
         ispinor=1
         ispinor1=1
       else
         write(message,'(2a)') " BUG in m_energy: nsploop should be equal to 1, 2 or 4"
         call wrtout(std_out,message,'COLL')
       endif
! Initialize noccmmp
       do im1 = 1 , ldim
         do im = 1 ,  ldim
            noccmmp(1,im,im1,idijeff)=real(green%occup%matlu(iatom)%mat(im1,im,isppol,ispinor,ispinor1))
            noccmmp(2,im,im1,idijeff)=aimag(green%occup%matlu(iatom)%mat(im1,im,isppol,ispinor,ispinor1))
!            noccmmp(1,im,im1,idijeff)=real(green%occup%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
!            noccmmp(2,im,im1,idijeff)=imag(green%occup%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
         enddo
       enddo
! Compute nocctot
       if(green%has_charge_matlu_solver/=2) then
         do im1=1,ldim
           if(nsploop==4) then
             if(idijeff<=2) then
               nocctot(1)=nocctot(1)+&
&                noccmmp(1,im1,im1,idijeff)
             endif
           else if (nsploop==2.or.nsploop==1) then
             nocctot(idijeff)=nocctot(idijeff)+&
&              noccmmp(1,im1,im1,idijeff)
           end if
         enddo
       else
         if(nsploop==4) then
           nocctot(1)=green%charge_matlu_solver(iatom,2) !  total nb of elec for nspinor=2 is (iatom,2) !!
           nocctot(2)=zero
           nocctot(3)=zero
           nocctot(4)=zero
          else if (nsploop==2) then
           nocctot(1)=green%charge_matlu_solver(iatom,1) !  first spin
           nocctot(2)=green%charge_matlu_solver(iatom,2) !  second one
          else if (nsploop==1) then
           nocctot(1)=green%charge_matlu_solver(iatom,1)
         end if
       endif
     enddo

     xe1=e_dc
     xe2=e_ee
    ! write(std_out,*)" nocctot(1)",nocctot(1),green%charge_matlu_solver(iatom,1)
     edftumdc = zero
     edftumdcdc = zero
     if ( present(renorm) ) then
       upawu = one
       jpawu = renorm(iatom)
       prt_pawuenergy=0
     else
       upawu = pawtab_%upawu
       jpawu = pawtab_%jpawu
       prt_pawuenergy=3
     end if

     call pawuenergy(iatom,edftumdc,edftumdcdc,noccmmp,nocctot,prt_pawuenergy,pawtab_,&
&                    dmft_dc=paw_dmft%dmft_dc,e_ee=e_ee,e_dc=e_dc,e_dcdc=e_dcdc,&
&                    u_dmft=upawu,j_dmft=jpawu)

     energies_dmft%e_dc(iatom)=e_dc-xe1
     energies_dmft%e_hu_dftu(iatom)=e_ee-xe2

     ABI_DEALLOCATE(noccmmp)
     ABI_DEALLOCATE(nocctot)
   endif ! lpawu/=-1
 enddo

! - gather results
! -----------------------------------------------------------------------
 energies_dmft%e_dc_tot=e_dc ! todo_ab: here or not ?
 energies_dmft%e_hu_dftu_tot=e_ee

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
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine compute_noninterentropy(cryst_struc,green,paw_dmft)

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
       eig=paw_dmft%eigen_dft(isppol,ikpt,ib)
       fermi=paw_dmft%fermie_dft
       fermi=paw_dmft%fermie
       occ1=occup_fd(eig,fermi,paw_dmft%temp)
       occ2=green%occup%ks(isppol,ikpt,ib,ib)
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

!!****f* m_energy/occup_fd
!! NAME
!! occup_fd
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

 function occup_fd(eig,fermie,temp)

!Arguments ------------------------------------
!type
! Integrate analytic tail 1/(iw-mu)
 real(dp),intent(in) :: eig,fermie,temp
 real(dp) :: occup_fd
!Local variables-------------------------------
! *********************************************************************

 if((eig-fermie) > zero) then
   occup_fd=exp(-(eig-fermie)/temp)/(one+exp(-(eig-fermie)/temp))
 else
   occup_fd=one/(one+exp((eig-fermie)/temp))
 endif

 end function occup_fd

END MODULE m_energy
!!***

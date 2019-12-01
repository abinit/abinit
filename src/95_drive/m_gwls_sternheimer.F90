!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gwls_sternheimer
!! NAME
!!   m_gwls_sternheimer
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2019 ABINIT group (JLJ, BR, MC)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_gwls_sternheimer

 implicit none

 private
!!***

 public :: gwls_sternheimer
!!***

contains
!!***

!!****f* ABINIT/gwls_sternheimer
!! NAME
!! gwls_sternheimer
!!
!! FUNCTION
!! .
!!
!! INPUTS
!! dtset <type(dataset_type)>=all input variables in this dataset
!!
!! OUTPUT
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      cleanupvalencewavefunctions,close_timing_log
!!      compute_correlations_no_model_shift_lanczos
!!      compute_correlations_shift_lanczos
!!      compute_exchange_and_correlation_energies,destroy_h
!!      driver_generateprintdielectriceigenvalues,gstateimg
!!      preparevalencewavefunctions,setup_timing_log,timab
!!
!! SOURCE

subroutine gwls_sternheimer(acell_img,amu_img,codvsn,cpui,dtfil,dtset,etotal_img,fcart_img,&
&                           fred_img,iexit,intgres_img,mixalch_img,mpi_enreg,nimage,npwtot,occ_img,&
&                           pawang,pawrad,pawtab,psps,rprim_img,strten_img,vel_cell_img,vel_img,xred_img,&
&                           filnam,filstat,idtset,jdtset,ndtset) ! optional arguments

use m_gwls_utility,                   only : master_debug, files_status_new, files_status_old
! use m_gwls_wf,                        only : norm_k
use m_gwls_hamiltonian,               only : destroy_H, exchange, g_to_r, eig
use m_gwls_TimingLog
use m_gwls_valenceWavefunctions
use m_gwls_ComputeCorrelationEnergy
use m_gwls_GenerateEpsilon
use m_dtset

use defs_basis
use defs_wvltypes
use m_pawang
use m_pawrad
use m_pawtab
use m_abicore
use m_errors
use m_dtfil

use defs_datatypes, only : pseudopotential_type
use defs_abitypes, only : MPI_type
use m_time,      only : timab
use m_gstateimg, only : gstateimg

!Arguments ------------------------------------
!scalars
integer,intent(in) :: nimage
integer,optional,intent(in) :: idtset,ndtset
integer,intent(inout) :: iexit
real(dp),intent(in) :: cpui
character(len=6),intent(in) :: codvsn
character(len=fnlen),optional,intent(in) :: filstat
type(MPI_type),intent(inout) :: mpi_enreg
type(datafiles_type),target,intent(inout) :: dtfil
type(dataset_type),intent(inout) :: dtset
type(pawang_type),intent(inout) :: pawang
type(pseudopotential_type),intent(inout) :: psps
!arrays
integer,optional,intent(in) :: jdtset(:)
integer,intent(out) :: npwtot(dtset%nkpt)
character(len=fnlen),optional,intent(in) :: filnam(:)
real(dp), intent(out) :: etotal_img(nimage),fcart_img(3,dtset%natom,nimage)
real(dp), intent(out) :: fred_img(3,dtset%natom,nimage),strten_img(6,nimage)
real(dp), intent(out) :: intgres_img(dtset%nspden,dtset%natom,nimage)
real(dp),intent(inout) :: acell_img(3,nimage),amu_img(dtset%ntypat,nimage),occ_img(dtset%mband*dtset%nkpt*dtset%nsppol,nimage)
real(dp),intent(inout) :: mixalch_img(dtset%npspalch,dtset%ntypalch,nimage)
real(dp),intent(inout) :: rprim_img(3,3,nimage),vel_cell_img(3,3,nimage),vel_img(3,dtset%natom,nimage)
real(dp),intent(inout) :: xred_img(3,dtset%natom,nimage)
type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
!Local variables ------------------------------

integer  :: e_index
real(dp) :: exchange_energy
real(dp) :: vxc_energy
real(dp) :: DFT_band_energy
 type(wvl_data) :: wvl

! timing
real(dp) :: tsec(2)
integer :: GWLS_TIMAB, OPTION_TIMAB
!END DEBUG
! *************************************************************************

! initialize to zero

! start new increment of the GWLS timer
 GWLS_TIMAB   = 1501
 OPTION_TIMAB = 1
 call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


 master_debug = .false.
!master_debug = .true.

!Governs the status of all opened files formerly with the status 'new' harcoded.
!It is initialized in the m_gwls_utility module as 'new'; this default is overriden here if desired.
 files_status_new = 'unknown'
!Governs the status of all opened files formerly with the status 'old' harcoded.
!It is initialized in the m_gwls_utility module as 'old'; this default is overriden here if desired.
 files_status_old = 'unknown'

! Test the input to make sure it is consistent with what
! the code can do
!call test_input(dtset2,psps2)

!Initializing variables and build the Hamiltonian using build_H

 GWLS_TIMAB   = 1521
 OPTION_TIMAB = 1
 call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


 call gstateimg(acell_img,amu_img,codvsn,cpui,dtfil,dtset,etotal_img,fcart_img,&
& fred_img,iexit,intgres_img,mixalch_img,mpi_enreg,nimage,npwtot,occ_img,&
& pawang,pawrad,pawtab,psps,rprim_img,strten_img,vel_cell_img,vel_img,wvl,xred_img,&
& filnam,filstat,idtset,jdtset,ndtset) ! optional arguments

 OPTION_TIMAB = 2
 call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


! tabulate the valence wavefunctions, to be used throughout the code!
! TODO : put in build_H and destroy_H
 GWLS_TIMAB   = 1522
 OPTION_TIMAB = 1
 call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

 call prepareValenceWavefunctions()

 OPTION_TIMAB = 2
 call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)



 e_index         = dtset%gwls_band_index
 DFT_band_energy = eig(e_index)

 write(std_out,10) '                               '
 write(std_out,10) ' GWLS RESULTS                  '
 write(std_out,10) ' ----------------------------- '
 write(std_out,12) ' For orbital |psi_e>, with e : ',e_index
 write(std_out,14) ' DFT eigenenergy             : ',DFT_band_energy,' Ha = ',DFT_band_energy*Ha_eV,' eV'
 flush(std_out)

 write(ab_out,10) '                               '
 write(ab_out,10) ' GWLS RESULTS                  '
 write(ab_out,10) ' ----------------------------- '
 write(ab_out,12) ' For orbital |psi_e>, with e : ',e_index
 write(ab_out,14) ' DFT eigenenergy             : ',DFT_band_energy,' Ha = ',DFT_band_energy*Ha_eV,' eV'

 if (dtset%gwls_exchange /= 0) then

   GWLS_TIMAB   = 1502
   OPTION_TIMAB = 1
   call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

   call compute_Exchange_and_Correlation_energies(e_index, exchange_energy, Vxc_energy)

   OPTION_TIMAB = 2
   call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)




   write(std_out,14) ' <psi_e |   V_xc  | psi_e>   : ',vxc_energy     ,' Ha = ',vxc_energy     *Ha_eV,' eV'
   write(std_out,14) ' <psi_e | Sigma_x | psi_e>   : ',exchange_energy,' Ha = ',exchange_energy*Ha_eV,' eV'
   flush(std_out)

   write(ab_out,14) ' <psi_e |   V_xc  | psi_e>   : ',vxc_energy     ,' Ha = ',vxc_energy     *Ha_eV,' eV'
   write(ab_out,14) ' <psi_e | Sigma_x | psi_e>   : ',exchange_energy,' Ha = ',exchange_energy*Ha_eV,' eV'
 else
   exchange_energy  = zero
   Vxc_energy       = zero

 end if

 if (dtset%gwls_correlation == 3) then

   call setup_timing_log()

   GWLS_TIMAB   = 1503
   OPTION_TIMAB = 1

   call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

   call compute_correlations_shift_lanczos(dtset, exchange_energy, Vxc_energy,master_debug)

   OPTION_TIMAB = 2
   call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


   call close_timing_log()

 else if (dtset%gwls_correlation == 4) then
   call setup_timing_log()
   call compute_correlations_no_model_shift_lanczos(dtset, exchange_energy, Vxc_energy,master_debug)
   call close_timing_log()

 else if (dtset%gwls_correlation == 5) then
   call setup_timing_log()
   call Driver_GeneratePrintDielectricEigenvalues(dtset)
   call close_timing_log()


 end if

 call cleanupValenceWavefunctions()

 call destroy_H()

 GWLS_TIMAB   = 1501
 OPTION_TIMAB = 2
 call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


 10 format(A)
 12 format(A,I6)
 14 format(A,ES24.16,A,F16.8,A)

end subroutine gwls_sternheimer
!!***

end module m_gwls_sternheimer
!!***

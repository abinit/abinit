!!****m* ABINIT/m_rttddft_driver
!! NAME
!!  m_rttddft_driver
!!
!! FUNCTION
!!  Real-Time Time Dependent DFT (RT-TDDFT)
!!
!! COPYRIGHT
!!  Copyright (C) 2021-2026 ABINIT group (FB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_rttddft_driver

 use defs_basis
 use defs_abitypes,        only: MPI_type
 use defs_datatypes,       only: pseudopotential_type
 use m_dtfil,              only: datafiles_type
 use m_dtset,              only: dataset_type
 use m_errors,             only: msg_hndl
 use m_pawang,             only: pawang_type
 use m_pawrad,             only: pawrad_type
 use m_pawtab,             only: pawtab_type
 use m_rttddft_output,     only: rttddft_output
 use m_rttddft_tdks,       only: tdks_type
 use m_rttddft_propagate,  only: rttddft_propagate_ele
 use m_rttddft_properties, only: rttddft_calc_density, &
                               & rttddft_calc_etot,    &
                               & rttddft_calc_current
 use m_specialmsg,         only: wrtout
 use m_time,               only: cwtime

 implicit none

 private
!!***

 public :: rttddft
!!***

contains
!!***

!!****f* m_rttddft_driver/rttddft
!! NAME
!!  rttddft
!!
!! FUNCTION
!!  Primary routine for real-time TDDFT calculations.
!!  1) Initialization: create main tdks (Time-Dependent Kohn-Sham) object
!!    - intialize various important parameters
!!    - read intial KS orbitals in WFK file
!!    - compute initial density from KS orbitals
!!  2) Propagation loop (rttddft_propagate):
!!    for i = 1, ntime
!!       - propagate KS orbitals
!!       - propagate nuclei if requested (Erhenfest dynamics)
!!       - compute new density
!!       - Compute and print requested properties
!!  3) Final printout and finalize
!!
!! INPUTS
!!  codvsn = code version
!!  dtfil <type datafiles_type> = infos about file names
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  pawang <type(pawang_type)> = paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)> = paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)> = paw tabulated starting data
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!
!! OUTPUT
!!
!! SOURCE
subroutine rttddft(codvsn,dtfil,dtset,mpi_enreg,pawang,pawrad,pawtab,psps)

 !Arguments ------------------------------------
 !scalars
 character(len=8),           intent(in)    :: codvsn
 type(dataset_type),         intent(inout) :: dtset
 type(datafiles_type),       intent(inout) :: dtfil
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pawang_type),          intent(inout) :: pawang
 type(pseudopotential_type), intent(inout) :: psps
 !arrays
 type(pawrad_type), target,  intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type), target,  intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

 !Local variables-------------------------------
 !scalars
 character(len=500)   :: msg
 integer              :: istep
 real(dp)             :: integrated_density, stability
 type(tdks_type)      :: tdks
 !arrays
 real(dp)             :: cpu, wall, gflops

! ***********************************************************************

 write(msg,'(7a)') ch10,'---------------------------------------------------------------------------',&
                 & ch10,'----------------   Starting real-time time dependent DFT   ----------------',&
                 & ch10,'---------------------------------------------------------------------------',ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)
 write(msg,'(7a)') ch10,' RT-TDDFT is under active development and should thus be used with caution ',&
                 & ch10,'---------------------------------------------------------------------------',ch10
 if (do_write_log) call wrtout(std_out,msg)

 !** 1) Initialization: create main tdks (Time-Dependent Kohn-Sham) object
 write(msg,'(3a)') ch10,'---------------------------   Initialization   ----------------------------',ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

 call tdks%init(codvsn,dtfil,dtset,mpi_enreg,pawang,pawrad,pawtab,psps)

 !Compute initial electronic density
 call rttddft_calc_density(dtset,mpi_enreg,psps,tdks)

 !Compute current at t=0 (or t-dt) and update vector potential accordingly
 if (dtset%prtcurrent/=0 .or. tdks%tdef%induced_vecpot) then
   call rttddft_calc_current(tdks,dtset,dtfil,psps,mpi_enreg)
 end if

 !** 2) Propagation loop
 write(msg,'(3a)') ch10,'-------------------------   Starting propagation   ------------------------',ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

 do istep = tdks%first_step, tdks%first_step+tdks%ntime-1

   call cwtime(cpu,wall,gflops,"start")

   !*** Perform electronic step ***
   !Compute new WF at time t and energy contribution at time t-dt
   call rttddft_propagate_ele(dtset,istep,mpi_enreg,psps,tdks)

   !Compute total energy at time t-dt
   call rttddft_calc_etot(dtset,tdks%energies,tdks%etot,tdks%occ)

   !Update electric field and vector potential value
   call tdks%tdef%update(dtset,mpi_enreg,istep*tdks%dt,tdks%rprimd,tdks%gprimd,tdks%kg, &
                       & psps%mpsang,tdks%npwarr,tdks%ylm,tdks%ylmgr,tdks%current)

   !Compute new electronic density at t
   call rttddft_calc_density(dtset,mpi_enreg,psps,tdks)

   !Compute current at t
   if (dtset%prtcurrent/=0 .or. tdks%tdef%induced_vecpot) then
      call rttddft_calc_current(tdks,dtset,dtfil,psps,mpi_enreg)
   end if

   !Compute and output useful electronic values
   call rttddft_output(dtfil,dtset,istep,mpi_enreg,psps,tdks)

   !Test of stability
   integrated_density = sum(tdks%rhor(:,1))*tdks%ucvol/tdks%nfftf
   stability = abs(integrated_density-dtset%nelect)/dtset%nelect
   if (stability > 0.1_dp) then
      write(msg,"(3a)") "The integrated density has changed by more than 10%!", ch10, &
                     &  "The integration is unstable, you should probably decrease the timestep dtele."
      ABI_ERROR(msg)
   else if (stability > 0.05_dp) then
      write(msg,"(3a)") "The integrated density has changed by more than 5%!", ch10, &
                     &  "You should be careful, the integration might be unstable"
      ABI_WARNING(msg)
   end if

   !TODO: *** Perform nuclear step ***
   ! For Ehrenfest dynamics
   !call rttddft_propagate_nuc(dtset,istep,mpi_enreg,psps,tdks)

   call cwtime(cpu,wall,gflops,"stop")
   write(msg,'(a,2f8.2,a)') '- Time - cpu, wall (sec):', cpu, wall, ch10
   call wrtout(ab_out,msg)
   if (do_write_log) call wrtout(std_out,msg)

 end do

 !** 3) Final Output and free memory
 write(msg,'(7a)') ch10,'---------------------------------------------------------------------------',&
                 & ch10,'----------------   Finished real-time time dependent DFT   ----------------',&
                 & ch10,'---------------------------------------------------------------------------',ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

 call tdks%free(dtset,mpi_enreg,psps)

end subroutine rttddft
!!***

end module m_rttddft_driver
!!***

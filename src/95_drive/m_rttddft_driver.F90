!!****m* ABINIT/m_rttddft_driver
!! NAME
!!  m_rttddft
!!
!! FUNCTION
!!  Real-time Time Dependent DFT (RT-TDDFT)
!!
!! COPYRIGHT
!!  Copyright (C) 2021 ABINIT group (FB, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
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

module m_rttddft_driver

 use defs_basis
 use defs_abitypes,       only: MPI_type
 use defs_datatypes,      only: pseudopotential_type

 use m_dtfil,             only: datafiles_type
 use m_dtset,             only: dataset_type
 use m_pawang,            only: pawang_type
 use m_pawrad,            only: pawrad_type
 use m_pawtab,            only: pawtab_type
 use m_rttddft,           only: rttddft_calc_density, &
                              & rttddft_calc_etot 
 use m_rttddft_output,    only: rttddft_output
 use m_rttddft_types,     only: tdks_type
 use m_rttddft_propagate, only: rttddft_propagate_ele
 use m_specialmsg,        only: wrtout

 implicit none

 private
!!***

 public :: rttddft
!!***

contains

!!****f* m_rttddft_driver/rttddft
!! NAME
!!  rttddft
!!
!! FUNCTION
!!  Primary routine that drives real-time TDDFT calculations.
!!  1) Initialization: create main tdks (Time-Dependent Kohn-Sham) object
!!    - intialize various important parameters
!!    - read KS orbitals in WFK file
!!    - compute initial density from KS orbitals
!!  2) Propagation loop (rttddft_propagate):
!!    for i = 1, nstep
!!       - propagate KS orbitals
!!       - propagate nuclei if requested (RTDDFT+Erhenfest dynamics)
!!       - compute new density and occupation
!!       - print requested quantities
!!  3) Final printout and finalize
!!
!! INPUTS
!!  codvsn = code version
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  dtfil <type datafiles_type> = infos about file names, file unit numbers
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  pawang <type(pawang_type)> = paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)> = paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)> = paw tabulated starting data
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  m_driver
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft(codvsn, dtfil, dtset, mpi_enreg, pawang, pawrad, pawtab, psps)

 implicit none

 !Arguments ------------------------------------
 !scalars
 character(len=8),           intent(in)    :: codvsn
 type(dataset_type),         intent(inout) :: dtset
 type(datafiles_type),       intent(inout) :: dtfil
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pawang_type),          intent(inout) :: pawang
 type(pseudopotential_type), intent(inout) :: psps
 !arrays
 type(pawrad_type),          intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),          intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

 !Local variables-------------------------------
 !scalars
 character(len=500)   :: msg
 integer              :: istep
 type(tdks_type)      :: tdks
 !arrays

! ***********************************************************************

 write(msg,'(7a)') ch10,'---------------------------------------------------------------------------',&
                 & ch10,'----------------   Starting real-time time dependent DFT   ----------------',&
                 & ch10,'---------------------------------------------------------------------------',ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

 !** 1) Initialization: create main tdks (Time-Dependent Kohn-Sham) object
 write(msg,'(3a)') ch10,'---------------------------   Initialization   ----------------------------',ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

 call tdks%init(codvsn,dtfil,dtset,mpi_enreg,pawang,pawrad,pawtab,psps)

 !Output useful values
 call rttddft_output(dtfil,dtset,0,mpi_enreg,tdks)

 !** 2) Propagation loop (rttddft_propagate):
 write(msg,'(3a)') ch10,'-------------------------   Starting propagation   ------------------------',ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)
 
 do istep = 1, tdks%ntime

   !FB TODO: If Ehrenfest perform nuclear step here
   !call rttddft_propagate_nuc(dtset,istep,mpi_enreg,psps,tdks)

   !Perform electronic step
   call rttddft_propagate_ele(dtset,istep,mpi_enreg,psps,tdks)

   !Calc total energy
   call rttddft_calc_etot(dtset,tdks%energies,tdks%etot)

   !Calc new electronic density 
   call rttddft_calc_density(dtset,mpi_enreg,psps,tdks)

   !Output useful values
   call rttddft_output(dtfil,dtset,istep,mpi_enreg,tdks)

 end do

 !** 3) Final Output and free memory
 write(msg,'(7a)') ch10,'---------------------------------------------------------------------------',&
                 & ch10,'----------------   Finished real-time time dependent DFT   ----------------',&
                 & ch10,'---------------------------------------------------------------------------',ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

 call tdks%free(dtset,mpi_enreg,psps)

end subroutine rttddft

end module m_rttddft_driver
!!***

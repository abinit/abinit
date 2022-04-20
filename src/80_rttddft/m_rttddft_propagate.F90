!!****m* ABINIT/m_rttddft_propagate
!! NAME
!!  m_rttddft_propagate
!!
!! FUNCTION
!!  Driver to perform electronic or nuclear step
!!
!! COPYRIGHT
!!  Copyright (C) 2021-2022 ABINIT group (FB)
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

module m_rttddft_propagate

 use defs_basis
 use defs_abitypes,         only: MPI_type
 use defs_datatypes,        only: pseudopotential_type
 
 use m_dtset,               only: dataset_type
 use m_errors,              only: msg_hndl
 use m_hamiltonian,         only: gs_hamiltonian_type
 use m_rttddft,             only: rttddft_setup_ele_step
 use m_rttddft_propagators, only: rttddft_propagator_er, &
                                & rttddft_propagator_emr
 use m_rttddft_tdks,        only: tdks_type
 use m_specialmsg,          only: wrtout

 implicit none

 private
!!***

 public :: rttddft_propagate_ele
 public :: rttddft_propagate_nuc
!!***

contains 
!!***

!!****f* m_rttddft/rttddft_propagate_ele
!!
!! NAME
!!  rttddft_propagate_ele
!!
!! FUNCTION
!!  Main subroutine to propagate time-dependent KS orbitals
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_propagate_ele(dtset, istep, mpi_enreg, psps, tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                    intent(in)    :: istep
 type(dataset_type),         intent(inout) :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 type(tdks_type),            intent(inout) :: tdks
 
 !Local variables-------------------------------
 !scalars
 character(len=500)        :: msg
 type(gs_hamiltonian_type) :: gs_hamk
 !arrays
 
! ***********************************************************************

 write(msg,'(a,a,i5)') ch10,'--- Iteration',istep
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

 ! Update various quantities after a nuclear step
 if (dtset%ionmov /= 0) call rttddft_setup_ele_step(dtset,mpi_enreg,psps,tdks)

 ! Propagate cg
 select case (dtset%td_propagator) 
   case(0)
      call rttddft_propagator_er(dtset,gs_hamk,istep,mpi_enreg,psps,tdks,calc_properties=.true.)
   case(1)
      call rttddft_propagator_emr(dtset,gs_hamk,istep,mpi_enreg,psps,tdks)  
   case default
      write(msg,"(a)") "Unknown Propagator - check the value of td_propagator"
      ABI_ERROR(msg)
 end select

end subroutine rttddft_propagate_ele
!!***

!!****f* m_rttddft/rttddft_propagate_nuc
!!
!! NAME
!!  rttddft_propagate_nuc
!!
!! FUNCTION
!!  Main subroutine to propagate nuclei using Ehrenfest dynamics
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_driver/rttddft
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_propagate_nuc(dtset, istep, mpi_enreg, psps, tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 type(tdks_type),           intent(inout) :: tdks
 integer,                    intent(in)    :: istep
 type(dataset_type),         intent(in)    :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 
 !Local variables-------------------------------
 !scalars
 character(len=500)   :: msg
 !arrays
 
! ***********************************************************************

 write(msg,'(2a,i5,a)') ch10,'--- Iteration',istep,ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

end subroutine rttddft_propagate_nuc
!!***

end module m_rttddft_propagate
!!***

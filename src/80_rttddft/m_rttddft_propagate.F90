!!****m* ABINIT/m_rttddft_propagate
!! NAME
!!  m_rttddft_propagate
!!
!! FUNCTION
!!  Contains a driver to propagate the KS orbitals 
!!  and the nuclei (Ehrenfest) if required
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
!!  m_rttddft_driver
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
 use m_rttddft_propagators, only: rttddft_propagator_er
 use m_rttddft_types,       only: tdks_type 
 use m_specialmsg,          only: wrtout
 use m_symtk,               only: symmetrize_xred

 implicit none

 private
!!***

 public :: rttddft_propagate_ele
 public :: rttddft_propagate_nuc
!!***

contains 

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
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
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
 character(len=500)             :: msg
 type(gs_hamiltonian_type)      :: gs_hamk
 !arrays
 
! ***********************************************************************

 write(msg,'(2a,i5,a)') ch10,'--- Iteration',istep,ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

 ! Init/Update various quantities before performing propagating KS orbitals
 call rttddft_setup_ele_step(dtset,gs_hamk,istep,mpi_enreg,psps,tdks)

 ! Propagate cg
 select case (dtset%td_propagator) 
   case(0)
      call rttddft_propagator_er(dtset,gs_hamk,istep,mpi_enreg,psps,tdks)
   case(1)
     ! call rttddft_propagator_emr(tdks,dtset,istep,mpi_enreg,psps)  
   case default
      write(msg,"(a,a)") "Unknown Propagator - check the value of td_propagator", ch10
      ABI_ERROR(msg)
 end select

 end subroutine rttddft_propagate_ele

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
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
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

 ! FB: Should we do this? 
 ! Eventually symmetrize atomic coordinates over space group elements:
 call symmetrize_xred(dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,tdks%xred,indsym=tdks%indsym)

 end subroutine rttddft_propagate_nuc

end module m_rttddft_propagate
!!***

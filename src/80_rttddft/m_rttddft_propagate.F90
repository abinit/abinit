!!****m* ABINIT/m_rttddft_propagate
!! NAME
!!  m_rttddft_propagate
!!
!! FUNCTION
!!  Contains various subroutines to propagate the KS 
!!  orbitals and potentially also the nuclei in RTTDDFT
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

module m_rttddft_propagate

 use defs_basis
 use defs_abitypes,     only: MPI_type
 use defs_datatypes,    only: pseudopotential_type
 
 use m_rttddft_types,   only: tdks_type 
 use m_specialmsg,      only: wrtout

 implicit none

 private
!!***

 public :: rttddft_propagate_ele

contains 

!!****f* m_rttddft/rttddft_propagate_ele
!!
!! NAME
!! rttddft_propagate_ele
!!
!! FUNCTION
!! Main subroutine to propagate time-dependent KS orbitals
!!
!! INPUTS
!! tdks <class(tdks_type)> = the tdks object to initialize
!! mpi_enreg <MPI_type> = MPI-parallelisation information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!! m_rttddft_driver
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_propagate_ele(tdks, itime, mpi_enreg, psps)

 implicit none

 !Arguments ------------------------------------
 !scalars
 class(tdks_type),           intent(inout) :: tdks
 integer,                    intent(in)    :: itime
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 
 !Local variables-------------------------------
 !scalars
 character(len=500)   :: msg
 !arrays
 
! ***********************************************************************

 write(msg,'(2a,i5,a)') ch10,'--- Iteration',itime,ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

 end subroutine rttddft_propagate_ele

end module m_rttddft_propagate
!!***

!!****m* ABINIT/m_rttddft_output
!! NAME
!!  m_rttddft_ouptut
!!
!! FUNCTION
!!  Manages most output of RT-TDDFT runs
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

module m_rttddft_output

 use defs_basis

 use m_errors,        only: msg_hndl, assert
 use m_rttddft_types, only: tdks_type
 use m_specialmsg,    only: wrtout
   
 implicit none

 private
!!***

 public :: rttddft_output
!!***

contains 

!!****f* m_rttddft_output/rttddft_output
!!
!! NAME
!!  rttddft_output
!!
!! FUNCTION
!!  Main output subroutine
!!
!! INPUTS
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
subroutine rttddft_output(tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 type(tdks_type),           intent(inout)        :: tdks
 !arrays
 
 !Local variables-------------------------------
 !scalars
 character(len=500)   :: msg
 !arrays

 write(msg,'(2a,f10.6)') ch10,'Integrated density (ie. total nb of electrons) = ', &
                    & SUM(tdks%rhor(:,1))*tdks%ucvol/tdks%nfftf
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

end subroutine rttddft_output

end module m_rttddft_output
!!***

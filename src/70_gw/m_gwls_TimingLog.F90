!!****m* ABINIT/m_gwls_TimingLog
!! NAME
!! m_gwls_TimingLog
!!
!! FUNCTION
!!  .
!!
!! COPYRIGHT
!! Copyright (C) 2009-2020 ABINIT group (JLJ, BR, MC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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


module m_gwls_TimingLog
!----------------------------------------------------------------------------------------------------
! This module will control the opening, writing and closing of a log file which will contain
! timing information. By making this a module, which can be called from anywhere, existing
! routines can easily write to this log without modifying their arguments.
!----------------------------------------------------------------------------------------------------
! local modules
use m_gwls_utility

! Abinit modules
use m_abicore
use defs_basis

use m_io_tools,  only : get_unit

implicit none
private
!!***

! Global quantities
integer, public  :: io_unit_timing_log

character(128)   :: timing_log_filename = "Timing_Log.dat"

logical :: head_node_timing
!!***

public :: close_timing_log, setup_timing_log
public :: write_block_lanczos_timing_log
public :: write_text_block_in_Timing_log
public :: write_timing_log
!!***

contains


!!****f* m_hamiltonian/setup_timing_log
!! NAME
!!  setup_timing_log
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_sternheimer
!!
!! CHILDREN
!!
!! SOURCE

subroutine setup_timing_log()
!--------------------------------------------------------------------------------
! This routine opens the timing log file.
!--------------------------------------------------------------------------------
use m_gwls_hamiltonian, only: mpi_enreg
implicit none


logical        :: file_exists
integer        :: i

! *************************************************************************


head_node_timing = (mpi_enreg%me == 0)


if (head_node_timing) then
  inquire(file=timing_log_filename,exist=file_exists)

  i = 0
  do while (file_exists)
  i = i+1
  write (timing_log_filename,'(A,I0,A)') "Timing_Log_",i,".dat"
  inquire(file=timing_log_filename,exist=file_exists)
  end do


  io_unit_timing_log = get_unit()
  open(io_unit_timing_log,file=timing_log_filename,status=files_status_new)
  write(io_unit_timing_log,10) ''
  write(io_unit_timing_log,10) '#============================================================================================'
  write(io_unit_timing_log,10) '# This file contains timing information for various routines from gw_sternheimer.'
  write(io_unit_timing_log,10) '# The purpose of this information is to establish which parts of the computation'
  write(io_unit_timing_log,10) '# are time consuming.'
  write(io_unit_timing_log,10) '#============================================================================================'
  write(io_unit_timing_log,10) '#   computation                                     time (seconds)       '
  write(io_unit_timing_log,10) '#============================================================================================'
  write(io_unit_timing_log,10) '#'
  flush(io_unit_timing_log)
end if

10 format(A)


end subroutine setup_timing_log
!!***

!!****f* m_hamiltonian/close_timing_log
!! NAME
!!  close_timing_log
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_sternheimer
!!
!! CHILDREN
!!
!! SOURCE

subroutine close_timing_log()
!--------------------------------------------------------------------------------
! This routine closes the timing log file.
!--------------------------------------------------------------------------------
implicit none
! *************************************************************************

if (head_node_timing) close(io_unit_timing_log)

end subroutine close_timing_log
!!***

!!****f* m_hamiltonian/write_text_block_in_Timing_log
!! NAME
!!  write_text_block_in_Timing_log
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_ComputeCorrelationEnergy,m_gwls_DielectricArray
!!
!! CHILDREN
!!
!! SOURCE

subroutine write_text_block_in_Timing_log(string)
!--------------------------------------------------------------------------------
! This routine opens the timing log file.
!--------------------------------------------------------------------------------
implicit none

character(256):: string
! *************************************************************************

if (head_node_timing) then
  write(io_unit_timing_log,10) trim(string)
  flush(io_unit_timing_log)
end if

10 format(A)

end subroutine write_text_block_in_Timing_log
!!***

!!****f* m_hamiltonian/write_timing_log
!! NAME
!!  write_timing_log
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_ComputeCorrelationEnergy,m_gwls_DielectricArray
!!      m_gwls_GenerateEpsilon
!!
!! CHILDREN
!!
!! SOURCE

subroutine write_timing_log(string,time)
!--------------------------------------------------------------------------------
! This routine opens the timing log file.
!--------------------------------------------------------------------------------
implicit none


character(256):: string
real(dp)     :: time
! *************************************************************************

if (head_node_timing) then
  write(io_unit_timing_log,20) trim(string),time
  flush(io_unit_timing_log)
end if


20 format(A,ES12.4)

end subroutine write_timing_log
!!***



!!****f* m_hamiltonian/write_block_lanczos_timing_log
!! NAME
!!  write_block_lanczos_timing_log
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_GWlanczos
!!
!! CHILDREN
!!
!! SOURCE

subroutine write_block_lanczos_timing_log(list_time,ntime)
!--------------------------------------------------------------------------------
! This routine writes the timing for the block lanczos routine.
!--------------------------------------------------------------------------------
use m_gwls_hamiltonian, only: mpi_enreg
implicit none


integer , intent(in) :: ntime
real(dp), intent(in) :: list_time(ntime)

logical        :: file_exists

integer, save  :: icounter = 0

integer        :: io_unit

character(128) :: block_lanczos_log_filename = "block_lanczos_timing.log"

logical        :: head_node

! *************************************************************************

head_node = (mpi_enreg%me == 0)
if (head_node) then

  inquire(file=block_lanczos_log_filename,exist=file_exists)

  io_unit = get_unit()

  if (.not. file_exists) then
    open(io_unit,file=block_lanczos_log_filename,status=files_status_new)
    write(io_unit,10) ''
    write(io_unit,10) '#==============================================================================================='
    write(io_unit,10) '# This file contains timing information for the block Lanczos algorithm.                        '
    write(io_unit,10) '# The purpose of this information is to establish which parts of the computation'
    write(io_unit,10) '# are time consuming.'
    write(io_unit,10) '#==============================================================================================='
    write(io_unit,10) '#   computation                                     time (seconds)       '
    write(io_unit,10) '#==============================================================================================='
    write(io_unit,10) '#'
  else 
    open(io_unit,file=block_lanczos_log_filename,position='append',status=files_status_old)
  end if


  icounter = icounter + 1 
  write(io_unit,10) ''
  write(io_unit,10) '#=========================='
  write(io_unit,15) '# Call number ', icounter 
  write(io_unit,10) '#=========================='
  write(io_unit,10) ''

  write(io_unit,20) '        Apply Matrix Function          :', list_time(1)
  write(io_unit,20) '        Compute Alpha                  :', list_time(2)
  write(io_unit,20) '        Update residual array (1-ZGEMM):', list_time(3)
  write(io_unit,20) '        Update residual array (2-ZGEMM):', list_time(4)
  write(io_unit,20) '        Orthogonalize                  :', list_time(5)
  write(io_unit,20) '        Extract QR factorization       :', list_time(6)
  write(io_unit,10) '-------------------------------------------------------'
  write(io_unit,20) '        Sum of the above               :', sum(list_time(1:6))
  write(io_unit,20) '        TOTAL ROUTINE TIME             :', list_time(7)


  close(io_unit)

end if 

10 format(A)
15 format(A,I5)
20 format(A,ES10.2)

end subroutine write_block_lanczos_timing_log
!!***

end module m_gwls_TimingLog
!!***

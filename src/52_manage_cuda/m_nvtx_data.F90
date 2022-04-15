!!****m* ABINIT/m_nvtx_data
!! NAME
!!  m_nvtx_data
!!
!! FUNCTION
!!  Small module to define data used to profile and trace abinit with nvtx library (Nvidia nsys).
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2000-2021 ABINIT group (MT)
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

module m_nvtx_data

#if defined(HAVE_GPU_CUDA) && defined(HAVE_GPU_NVTX_V3)
  use m_nvtx, only : nvtxStartRange, nvtxEndRange
#endif

  implicit none

  logical :: nvtx_activated = .false.

  integer, parameter :: NUMBER_OF_NVTX_REGIONS = 8
  character(len=32), dimension(NUMBER_OF_NVTX_REGIONS) :: nvtx_names
  integer          , dimension(NUMBER_OF_NVTX_REGIONS) :: nvtx_ids

  integer, parameter :: NVTX_MAIN_COMPUTATION = 1
  integer, parameter :: NVTX_SCF = 2
  integer, parameter :: NVTX_VTORHO = 3
  integer, parameter :: NVTX_LOBPCG1 = 4
  integer, parameter :: NVTX_LOBPCG2 = 5
  integer, parameter :: NVTX_CHEBFI1 = 6
  integer, parameter :: NVTX_CHEBFI2 = 7
  integer, parameter :: NVTX_ORTHO_WF = 8

contains

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine nvtx_init(activate)

    implicit none

    ! dummy variables
    logical :: activate

    nvtx_activated = activate

    nvtx_names = [character(len=32) :: "MAIN_COMPUTATION", "SCF", "VTORHO", &
         & "LOBPCG1", "LOBPCG2", &
         & "CHEBFI1", "CHEBFI2", &
         & "ORTHO_WF"]

    nvtx_ids(1) = NVTX_MAIN_COMPUTATION
    nvtx_ids(2) = NVTX_SCF
    nvtx_ids(3) = NVTX_VTORHO
    nvtx_ids(4) = NVTX_LOBPCG1
    nvtx_ids(5) = NVTX_LOBPCG2
    nvtx_ids(6) = NVTX_CHEBFI1
    nvtx_ids(7) = NVTX_CHEBFI2
    nvtX_ids(8) = NVTX_ORTHO_WF

  end subroutine nvtx_init

#if defined(HAVE_GPU_CUDA) && defined(HAVE_GPU_NVTX_V3)
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine abi_nvtx_start_range(id)

    implicit none

    ! dummy variables
    integer :: id

    if (nvtx_activated) then
       if (id .le. NUMBER_OF_NVTX_REGIONS) then
          call nvtxStartRange(nvtx_names(id),id)
       end if
    end if

  end subroutine abi_nvtx_start_range

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine abi_nvtx_end_range()

    implicit none

    if (nvtx_activated) then
       call nvtxEndRange()
    end if

  end subroutine abi_nvtx_end_range
#endif

end module m_nvtx_data

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

  integer, parameter :: NUMBER_OF_NVTX_REGIONS = 27
  character(len=32), dimension(NUMBER_OF_NVTX_REGIONS) :: nvtx_names
  integer          , dimension(NUMBER_OF_NVTX_REGIONS) :: nvtx_ids

  integer, parameter :: NVTX_MAIN_COMPUTATION = 1
  integer, parameter :: NVTX_SCF = 2
  integer, parameter :: NVTX_VTORHO = 3
  integer, parameter :: NVTX_VTOWFK = 4
  integer, parameter :: NVTX_LOBPCG1 = 5
  integer, parameter :: NVTX_LOBPCG2 = 6
  integer, parameter :: NVTX_CHEBFI1 = 7
  integer, parameter :: NVTX_CHEBFI2 = 8
  integer, parameter :: NVTX_ORTHO_WF = 9
  integer, parameter :: NVTX_GETGHC = 10
  integer, parameter :: NVTX_CHEBFI2_RR = 11 ! Rayleigh Ritz
  integer, parameter :: NVTX_CHEBFI2_RRQ = 12 ! Rayleigh Ritz Quotient
  integer, parameter :: NVTX_CHEBFI2_CORE = 13 ! core
  integer, parameter :: NVTX_CHEBFI2_NONLOP = 14
  integer, parameter :: NVTX_CTOCPRJ = 15
  integer, parameter :: NVTX_SCF_FOURWF = 16
  integer, parameter :: NVTX_MKRHO = 17
  integer, parameter :: NVTX_INVOVL = 18
  integer, parameter :: NVTX_INVOVL_NONLOP1 = 19
  integer, parameter :: NVTX_INVOVL_NONLOP2 = 20
  integer, parameter :: NVTX_INVOVL_INNER = 21
  integer, parameter :: NVTX_INVOVL_INNER_APPLY_BLOCK = 22
  integer, parameter :: NVTX_INVOVL_INNER_GEMM = 23
  integer, parameter :: NVTX_SUB_SPC_DIAGO = 24
  integer, parameter :: NVTX_CHEBFI2_NEXT_ORDER = 25
  integer, parameter :: NVTX_CHEBFI2_SWAP_BUF = 26
  integer, parameter :: NVTX_CHEBFI2_GET_AX_BX = 27

contains

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine nvtx_init(activate)

    implicit none

    ! dummy variables
    logical :: activate

    nvtx_activated = activate

    nvtx_names = [character(len=32) :: &
         & "MAIN_COMPUTATION", &
         & "SCF", &
         & "VTORHO", &
         & "VTOWFK", &
         & "LOBPCG1", &
         & "LOBPCG2", &
         & "CHEBFI1", &
         & "CHEBFI2", &
         & "ORTHO_WF", &
         & "GETGHC", &
         & "CHEBFI2_RR", &
         & "CHEBFI2_RRQ", &
         & "CHEBFI2_CORE", &
         & "CHEBFI2_NONLOP", &
         & "CTOCPRJ", &
         & "SCF_FOURWF", &
         & "MKRHO", &
         & "INVOVL", &
         & "INVOVL_NONLOP1", &
         & "INVOVL_NONLOP2", &
         & "INVOVL_INNER", &
         & "INVOVL_INNER_APPLY_BLOCK", &
         & "INVOVL_INNER_GEMM", &
         & "SUB_SPC_DIAGO", &
         & "CHEBFI2_NEXT_ORDER", &
         & "CHEBFI2_SWAP_BUF", &
         & "CHEBFI2_GET_AX_BX" &
         ]

    nvtx_ids(1) = NVTX_MAIN_COMPUTATION
    nvtx_ids(2) = NVTX_SCF
    nvtx_ids(3) = NVTX_VTORHO
    nvtx_ids(4) = NVTX_VTOWFK
    nvtx_ids(5) = NVTX_LOBPCG1
    nvtx_ids(6) = NVTX_LOBPCG2
    nvtx_ids(7) = NVTX_CHEBFI1
    nvtx_ids(8) = NVTX_CHEBFI2
    nvtx_ids(9) = NVTX_ORTHO_WF
    nvtx_ids(10) = NVTX_GETGHC
    nvtx_ids(11)= NVTX_CHEBFI2_RR
    nvtx_ids(12)= NVTX_CHEBFI2_RRQ
    nvtx_ids(13)= NVTX_CHEBFI2_CORE
    nvtx_ids(14)= NVTX_CHEBFI2_NONLOP
    nvtx_ids(15)= NVTX_CTOCPRJ
    nvtx_ids(16)= NVTX_SCF_FOURWF
    nvtx_ids(17)= NVTX_MKRHO
    nvtx_ids(18)= NVTX_INVOVL
    nvtx_ids(19)= NVTX_INVOVL_NONLOP1
    nvtx_ids(20)= NVTX_INVOVL_NONLOP2
    nvtx_ids(21)= NVTX_INVOVL_INNER
    nvtx_ids(22)= NVTX_INVOVL_INNER_APPLY_BLOCK
    nvtx_ids(23)= NVTX_INVOVL_INNER_GEMM
    nvtx_ids(24)= NVTX_SUB_SPC_DIAGO
    nvtx_ids(25)= NVTX_CHEBFI2_NEXT_ORDER
    nvtx_ids(26)= NVTX_CHEBFI2_SWAP_BUF
    nvtx_ids(27)= NVTX_CHEBFI2_GET_AX_BX

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
!!***

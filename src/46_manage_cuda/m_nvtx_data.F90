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

  integer, parameter :: NUMBER_OF_NVTX_REGIONS = 46
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
  integer, parameter :: NVTX_GETGHC_LOCPOT = 11
  integer, parameter :: NVTX_GETGHC_NLOCPOT = 12
  integer, parameter :: NVTX_GETGHC_KIN = 13
  integer, parameter :: NVTX_CHEBFI2_RR = 14 ! Rayleigh Ritz
  integer, parameter :: NVTX_CHEBFI2_RRQ = 15 ! Rayleigh Ritz Quotient
  integer, parameter :: NVTX_CHEBFI2_CORE = 16 ! core
  integer, parameter :: NVTX_CHEBFI2_NONLOP = 17
  integer, parameter :: NVTX_CTOCPRJ = 18
  integer, parameter :: NVTX_SCF_FOURWF = 19
  integer, parameter :: NVTX_MKRHO = 20
  integer, parameter :: NVTX_INVOVL = 21
  integer, parameter :: NVTX_INVOVL_PREP = 22
  integer, parameter :: NVTX_INVOVL_NONLOP1 = 23
  integer, parameter :: NVTX_INVOVL_NONLOP2 = 24
  integer, parameter :: NVTX_INVOVL_INNER = 25
  integer, parameter :: NVTX_INVOVL_POST1 = 26
  integer, parameter :: NVTX_INVOVL_POST2 = 27
  integer, parameter :: NVTX_INVOVL_POST3 = 28
  integer, parameter :: NVTX_SUB_SPC_DIAGO = 29
  integer, parameter :: NVTX_CHEBFI2_NEXT_ORDER = 30
  integer, parameter :: NVTX_CHEBFI2_SWAP_BUF = 31
  integer, parameter :: NVTX_CHEBFI2_GET_AX_BX = 32
  integer, parameter :: NVTX_VTOWFK_EXTRA1 = 33
  integer, parameter :: NVTX_VTOWFK_EXTRA2 = 34
  integer, parameter :: NVTX_VTORHO_EXTRA = 35
  integer, parameter :: NVTX_SCFCV_PAWKNHAT = 36
  integer, parameter :: NVTX_SCFCV_RHOTOV = 37
  integer, parameter :: NVTX_SCFCV_NEWRHO = 38
  integer, parameter :: NVTX_SCFCV_PAWDENPOT = 39
  integer, parameter :: NVTX_SCFCV_ETOTFOR = 40
  integer, parameter :: NVTX_SCFCV_SCPRQT = 41
  integer, parameter :: NVTX_SCFCV_DIJ = 42
  integer, parameter :: NVTX_SCFCV_SETVTR = 43
  integer, parameter :: NVTX_CHEBFI2_SQRT2 = 44
  integer, parameter :: NVTX_CHEBFI2_INIT = 45
  integer, parameter :: NVTX_INIT_INWFFIL = 46

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
         & "LOCPOT", &
         & "NLOCPOT", &
         & "KINETIC", &
         & "CHEBFI2_RR", &
         & "CHEBFI2_RRQ", &
         & "CHEBFI2_CORE", &
         & "CHEBFI2_NONLOP", &
         & "CTOCPRJ", &
         & "SCF_FOURWF", &
         & "MKRHO", &
         & "INVOVL", &
         & "INVOVL_PREP", &
         & "INVOVL_NONLOP1", &
         & "INVOVL_NONLOP2", &
         & "INVOVL_INNER", &
         & "INVOVL_POST1", &
         & "INVOVL_POST2", &
         & "INVOVL_POST3", &
         & "SUB_SPC_DIAGO", &
         & "CHEBFI2_NEXT_ORDER", &
         & "CHEBFI2_SWAP_BUF", &
         & "CHEBFI2_GET_AX_BX", &
         & "VTOWFK_EXTRA1", &
         & "VTOWFK_EXTRA2", &
         & "VTORHO_EXTRA", &
         & "SCFCV_PAWKNHAT", &
         & "SCFCV_RHOTOV", &
         & "SCFCV_NEWRHO", &
         & "SCFCV_PAWDENPOT", &
         & "SCFCV_ETOTFOR", &
         & "SCFCV_SCPRQT", &
         & "SCFCV_DIJ", &
         & "SCFCV_SETVTR", &
         & "CHEBFI2_SQRT2", &
         & "CHEBFI2_INIT", &
         & "INIT_INWFFIL" &
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
    nvtx_ids(10)= NVTX_GETGHC
    nvtx_ids(11)= NVTX_GETGHC_LOCPOT
    nvtx_ids(12)= NVTX_GETGHC_NLOCPOT
    nvtx_ids(13)= NVTX_GETGHC_KIN
    nvtx_ids(14)= NVTX_CHEBFI2_RR
    nvtx_ids(15)= NVTX_CHEBFI2_RRQ
    nvtx_ids(16)= NVTX_CHEBFI2_CORE
    nvtx_ids(17)= NVTX_CHEBFI2_NONLOP
    nvtx_ids(18)= NVTX_CTOCPRJ
    nvtx_ids(19)= NVTX_SCF_FOURWF
    nvtx_ids(20)= NVTX_MKRHO
    nvtx_ids(21)= NVTX_INVOVL
    nvtx_ids(22)= NVTX_INVOVL_PREP
    nvtx_ids(23)= NVTX_INVOVL_NONLOP1
    nvtx_ids(24)= NVTX_INVOVL_NONLOP2
    nvtx_ids(25)= NVTX_INVOVL_INNER
    nvtx_ids(26)= NVTX_INVOVL_POST1
    nvtx_ids(27)= NVTX_INVOVL_POST2
    nvtx_ids(28)= NVTX_INVOVL_POST3
    nvtx_ids(29)= NVTX_SUB_SPC_DIAGO
    nvtx_ids(30)= NVTX_CHEBFI2_NEXT_ORDER
    nvtx_ids(31)= NVTX_CHEBFI2_SWAP_BUF
    nvtx_ids(32)= NVTX_CHEBFI2_GET_AX_BX
    nvtx_ids(33)= NVTX_VTOWFK_EXTRA1
    nvtx_ids(34)= NVTX_VTOWFK_EXTRA2
    nvtx_ids(35)= NVTX_VTORHO_EXTRA
    nvtx_ids(36)= NVTX_SCFCV_PAWKNHAT
    nvtx_ids(37)= NVTX_SCFCV_RHOTOV
    nvtx_ids(38)= NVTX_SCFCV_NEWRHO
    nvtx_ids(39)= NVTX_SCFCV_PAWDENPOT
    nvtx_ids(40)= NVTX_SCFCV_ETOTFOR
    nvtx_ids(41)= NVTX_SCFCV_SCPRQT
    nvtx_ids(42)= NVTX_SCFCV_DIJ
    nvtx_ids(43)= NVTX_SCFCV_SETVTR
    nvtx_ids(44)= NVTX_CHEBFI2_SQRT2
    nvtx_ids(45)= NVTX_CHEBFI2_INIT
    nvtx_ids(46)= NVTX_INIT_INWFFIL

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

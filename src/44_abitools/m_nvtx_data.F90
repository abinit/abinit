!!****m* ABINIT/m_nvtx_data
!! NAME
!!  m_nvtx_data
!!
!! FUNCTION
!!  Small module to define data used to profile and trace abinit with nvtx library (Nvidia nsys).
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2000-2025 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_nvtx_data

#if defined(HAVE_GPU_MARKERS)
  use m_nvtx, only : nvtxStartRange, nvtxEndRange, nvtxProfilerStart, nvtxProfilerStop
#endif

  implicit none

  integer, parameter :: NUMBER_OF_NVTX_REGIONS = 117
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
  integer, parameter :: NVTX_VTOWFK_FOURWF = 19
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
  integer, parameter :: NVTX_INIT_INWFFIL2 = 47
  integer, parameter :: NVTX_CHEBFI2_RR_SCALE = 48 ! not used anymore
  integer, parameter :: NVTX_RR_HEGV = 49
  integer, parameter :: NVTX_RR_GEMM_1 = 50
  integer, parameter :: NVTX_RR_GEMM_2 = 51
  integer, parameter :: NVTX_INVOVL_INNER_GEMM = 52
  integer, parameter :: NVTX_CHEBFI2_TRANSPOSE = 53
  integer, parameter :: NVTX_CHEBFI2_GET_BM1X = 54
  integer, parameter :: NVTX_LOBPCG2_BLOCK = 55
  integer, parameter :: NVTX_LOBPCG2_LINE = 56
  integer, parameter :: NVTX_LOBPCG2_ORTHO_X_WRT = 57
  integer, parameter :: NVTX_LOBPCG2_RESIDUE = 58
  integer, parameter :: NVTX_LOBPCG2_RR = 59
  integer, parameter :: NVTX_B_ORTHO = 60
  integer, parameter :: NVTX_LOBPCG2_GET_AX_BX = 61
  integer, parameter :: NVTX_RR_HEEV = 62
  integer, parameter :: NVTX_FORSTRNPS = 63
  integer, parameter :: NVTX_FORSTR_NONLOP = 64
  integer, parameter :: NVTX_VTOWFK_NONLOP = 65
  integer, parameter :: NVTX_DFPT_SCF = 66
  integer, parameter :: NVTX_DFPT_VTORHO = 67
  integer, parameter :: NVTX_DFPT_VTOWFK = 68
  integer, parameter :: NVTX_DFPT_CGWF = 69
  integer, parameter :: NVTX_DFPT_CGWF_CORE = 70
  integer, parameter :: NVTX_DFPT_LOOP = 71
  integer, parameter :: NVTX_DFPT_NSTPAW = 72
  integer, parameter :: NVTX_DFPT_NSTDY = 73
  integer, parameter :: NVTX_DFPT_NSTWF = 74
  integer, parameter :: NVTX_DFPT_NSTWF_BAND = 75
  integer, parameter :: NVTX_GETGH1C = 76
  integer, parameter :: NVTX_GETGH2C = 77
  integer, parameter :: NVTX_DFPT_LOOP_DDK = 78
  integer, parameter :: NVTX_DFPT_LOOP_EFELD = 79
  integer, parameter :: NVTX_DFPT_LOOP_STRAIN = 80
  integer, parameter :: NVTX_DFPT_LOOP_PHONON = 81
  integer, parameter :: NVTX_FOURWF = 82
  integer, parameter :: NVTX_NONLOP = 83
  integer, parameter :: NVTX_DFPT_NSELT = 84
  integer, parameter :: NVTX_RESPFN = 85
  integer, parameter :: NVTX_DFPT_ELT = 86
  integer, parameter :: NVTX_DFPT_ATM2FFT = 87
  integer, parameter :: NVTX_DFPT_DYXC = 88
  integer, parameter :: NVTX_D2FRNL = 89
  integer, parameter :: NVTX_D2FRNL_KPT = 90
  integer, parameter :: NVTX_DFPT_RHOFERMI = 91
  integer, parameter :: NVTX_DFPT_WFKFERMI = 92
  integer, parameter :: NVTX_GETGSC = 93
  integer, parameter :: NVTX_DFPT_ACCRHO = 94
  integer, parameter :: NVTX_DFPT_MKRHO = 95
  integer, parameter :: NVTX_MAKE_INVOVL = 96
  integer, parameter :: NVTX_FORSTR = 97
  integer, parameter :: NVTX_FORCES = 98
  integer, parameter :: NVTX_STRESS = 99
  integer, parameter :: NVTX_DMFT_SOLVE = 100
  integer, parameter :: NVTX_DMFT_SOLVE_LOOP = 101
  integer, parameter :: NVTX_DMFT_IMPURITY_SOLVE = 102
  integer, parameter :: NVTX_DMFT_HUBBARD_ONE = 103
  integer, parameter :: NVTX_DMFT_DOWNFOLD_OPER = 104
  integer, parameter :: NVTX_DMFT_UPFOLD_OPER = 105
  integer, parameter :: NVTX_DMFT_INVERSE_OPER = 106
  integer, parameter :: NVTX_DMFT_COMPUTE_GREEN = 107
  integer, parameter :: NVTX_DMFT_COMPUTE_GREEN_BATCHED = 108
  integer, parameter :: NVTX_DMFT_COMPUTE_GREEN_LOOP = 109
  integer, parameter :: NVTX_DMFT_INTEGRATE_GREEN = 110
  integer, parameter :: NVTX_DMFT_FERMI_GREEN = 111
  integer, parameter :: NVTX_DMFT_COMPUTE_NB_ELEC = 112
  integer, parameter :: NVTX_DMFT_ADD_INT_FCT = 113
  integer, parameter :: NVTX_DMFT_SYM_MATLU = 114
  integer, parameter :: NVTX_DMFT_RW_SELF = 115
  integer, parameter :: NVTX_DMFT_SAVEOCC = 116
  integer, parameter :: NVTX_TRANSPOSER_MPI_ALL2ALL = 117

contains

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine nvtx_init()

    implicit none

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
         & "RAYLRITZ", &
         & "RAYLRITZ_Q", &
         & "CHEBFI2_CORE", &
         & "CHEBFI2_NONLOP", &
         & "CTOCPRJ", &
         & "VTOWFK_FOURWF", &
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
         & "NEXT_ORDER", &
         & "SWAP_BUF", &
         & "GET_AX_BX", &
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
         & "INIT_INWFFIL", &
         & "INIT_INWFFIL2", &
         & "RR_SCALE", &
         & "RR_HEGV", &
         & "RR_XNP", &
         & "RR_GEMM", &
         & "INVOVL_INNER_GEMM", &
         & "TRANSPOSE", &
         & "GET_BM1X", &
         & "LOBPCG2_BLOCK", &
         & "LOBPCG2_LINE", &
         & "LOBPCG2_ORTHO_X_WRT", &
         & "LOBPCG2_RESIDUE", &
         & "LOBPCG2_RR", &
         & "B_ORTHO", &
         & "GET_AX_BX", &
         & "RR_HEEV", &
         & "FORSTRNPS", &
         & "FORSTR_NONLOP", &
         & "VTOWFK_NONLOP", &
         & "DFPT_SCF", &
         & "DFPT_VTORHO", &
         & "DFPT_VTOWFK", &
         & "DFPT_CGWF", &
         & "DFPT_CGWF_CORE", &
         & "DFPT_LOOP", &
         & "DFPT_NSTPAW", &
         & "DFPT_NSTDY", &
         & "DFPT_NSTWF", &
         & "DFPT_NSTWF_BAND", &
         & "GETGH1C", &
         & "GETGH2C", &
         & "DFPT_LOOP_DDK", &
         & "DFPT_LOOP_EFELD", &
         & "DFPT_LOOP_STRAIN", &
         & "DFPT_LOOP_PHONON", &
         & "FOURWF", &
         & "NONLOP", &
         & "DFPT_NSELT", &
         & "RESPFN", &
         & "DFPT_ELT", &
         & "DFPT_ATM2FFT", &
         & "DFPT_DYXC", &
         & "D2FRNL", &
         & "D2FRNL_KPT", &
         & "DFPT_RHOFERMI", &
         & "DFPT_WFKFERMI", &
         & "GETGSC", &
         & "DFPT_ACCRHO", &
         & "DFPT_MKRHO", &
         & "MAKE_INVOVL", &
         & "FORSTR", &
         & "FORCES", &
         & "STRESS", &
         & "DMFT_SOLVE", &
         & "DMFT_SOLVE_LOOP", &
         & "DMFT_IMPURITY_SOLVE", &
         & "DMFT_HUBBARD_ONE", &
         & "DMFT_DOWNFOLD_OPER", &
         & "DMFT_UPFOLD_OPER", &
         & "DMFT_INVERSE_OPER", &
         & "DMFT_COMPUTE_GREEN", &
         & "DMFT_COMPUTE_GREEN_BATCHED", &
         & "DMFT_COMPUTE_GREEN_LOOP", &
         & "DMFT_INTEGRATE_GREEN", &
         & "DMFT_FERMI_GREEN", &
         & "DMFT_COMPUTE_NB_ELEC", &
         & "DMFT_ADD_INT_FCT", &
         & "DMFT_SYM_MATLU", &
         & "DMFT_RW_SELF", &
         & "DMFT_SAVEOCC", &
         & "TRANSPOSER_MPI_ALL2ALL" &
         & ]

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
    nvtx_ids(19)= NVTX_VTOWFK_FOURWF
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
    nvtx_ids(47)= NVTX_INIT_INWFFIL2
    nvtx_ids(48)= NVTX_CHEBFI2_RR_SCALE
    nvtx_ids(49)= NVTX_RR_HEGV
    nvtx_ids(50)= NVTX_RR_GEMM_1
    nvtx_ids(51)= NVTX_RR_GEMM_2
    nvtx_ids(52)= NVTX_INVOVL_INNER_GEMM
    nvtx_ids(53)= NVTX_CHEBFI2_TRANSPOSE
    nvtx_ids(54)= NVTX_CHEBFI2_GET_BM1X
    nvtx_ids(55)= NVTX_LOBPCG2_BLOCK
    nvtx_ids(56)= NVTX_LOBPCG2_LINE
    nvtx_ids(57)= NVTX_LOBPCG2_ORTHO_X_WRT
    nvtx_ids(58)= NVTX_LOBPCG2_RESIDUE
    nvtx_ids(59)= NVTX_LOBPCG2_RR
    nvtx_ids(60)= NVTX_B_ORTHO
    nvtx_ids(61)= NVTX_LOBPCG2_GET_AX_BX
    nvtx_ids(62)= NVTX_RR_HEEV
    nvtx_ids(63)= NVTX_FORSTRNPS
    nvtx_ids(64)= NVTX_FORSTR_NONLOP
    nvtx_ids(65)= NVTX_VTOWFK_NONLOP
    nvtx_ids(66)= NVTX_DFPT_SCF
    nvtx_ids(67)= NVTX_DFPT_VTORHO
    nvtx_ids(68)= NVTX_DFPT_VTOWFK
    nvtx_ids(69)= NVTX_DFPT_CGWF
    nvtx_ids(70)= NVTX_DFPT_CGWF_CORE
    nvtx_ids(71)= NVTX_DFPT_LOOP
    nvtx_ids(72)= NVTX_DFPT_NSTPAW
    nvtx_ids(73)= NVTX_DFPT_NSTDY
    nvtx_ids(74)= NVTX_DFPT_NSTWF
    nvtx_ids(75)= NVTX_DFPT_NSTWF_BAND
    nvtx_ids(76)= NVTX_GETGH1C
    nvtx_ids(77)= NVTX_GETGH2C
    nvtx_ids(78)= NVTX_DFPT_LOOP_DDK
    nvtx_ids(79)= NVTX_DFPT_LOOP_EFELD
    nvtx_ids(80)= NVTX_DFPT_LOOP_STRAIN
    nvtx_ids(81)= NVTX_DFPT_LOOP_PHONON
    nvtx_ids(82)= NVTX_FOURWF
    nvtx_ids(83)= NVTX_NONLOP
    nvtx_ids(84)= NVTX_DFPT_NSELT
    nvtx_ids(85)= NVTX_RESPFN
    nvtx_ids(86)= NVTX_DFPT_ELT
    nvtx_ids(87)= NVTX_DFPT_ATM2FFT
    nvtx_ids(88)= NVTX_DFPT_DYXC
    nvtx_ids(89)= NVTX_D2FRNL
    nvtx_ids(90)= NVTX_D2FRNL_KPT
    nvtx_ids(91)= NVTX_DFPT_RHOFERMI
    nvtx_ids(92)= NVTX_DFPT_WFKFERMI
    nvtx_ids(93)= NVTX_GETGSC
    nvtx_ids(94)= NVTX_DFPT_ACCRHO
    nvtx_ids(95)= NVTX_DFPT_MKRHO
    nvtx_ids(96)= NVTX_MAKE_INVOVL
    nvtx_ids(97)= NVTX_FORSTR
    nvtx_ids(98)= NVTX_FORCES
    nvtx_ids(99)= NVTX_STRESS
    nvtx_ids(100)=NVTX_DMFT_SOLVE
    nvtx_ids(101)=NVTX_DMFT_SOLVE_LOOP
    nvtx_ids(102)=NVTX_DMFT_IMPURITY_SOLVE
    nvtx_ids(103)=NVTX_DMFT_HUBBARD_ONE
    nvtx_ids(104)=NVTX_DMFT_DOWNFOLD_OPER
    nvtx_ids(105)=NVTX_DMFT_UPFOLD_OPER
    nvtx_ids(106)=NVTX_DMFT_INVERSE_OPER
    nvtx_ids(107)=NVTX_DMFT_COMPUTE_GREEN
    nvtx_ids(108)=NVTX_DMFT_COMPUTE_GREEN_BATCHED
    nvtx_ids(109)=NVTX_DMFT_COMPUTE_GREEN_LOOP
    nvtx_ids(110)=NVTX_DMFT_INTEGRATE_GREEN
    nvtx_ids(111)=NVTX_DMFT_FERMI_GREEN
    nvtx_ids(112)=NVTX_DMFT_COMPUTE_NB_ELEC
    nvtx_ids(113)=NVTX_DMFT_ADD_INT_FCT
    nvtx_ids(114)=NVTX_DMFT_SYM_MATLU
    nvtx_ids(115)=NVTX_DMFT_RW_SELF
    nvtx_ids(116)=NVTX_DMFT_SAVEOCC
    nvtx_ids(117)=NVTX_TRANSPOSER_MPI_ALL2ALL

  end subroutine nvtx_init

#if defined(HAVE_GPU_MARKERS)
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine abi_nvtx_start_range(id)

    implicit none

    ! dummy variables
    integer :: id

    if (id .le. NUMBER_OF_NVTX_REGIONS) then
       call nvtxStartRange(nvtx_names(id),id)
    end if

  end subroutine abi_nvtx_start_range

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine abi_nvtx_end_range()

    implicit none

    call nvtxEndRange()

  end subroutine abi_nvtx_end_range

#endif

end module m_nvtx_data
!!***

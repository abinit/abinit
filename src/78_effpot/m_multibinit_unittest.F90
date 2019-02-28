!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_multibinit_unittest
!! NAME
!!m_multibinit_unittest
!!
!! FUNCTION
!! This module contains subroutines for unit test
!! Instead of running multibinit, it will run unit tests.
!!
!! Datatypes:
!!
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
!! mb_test_main: run all unit test subroutines
!! mb_test1: 
!!
!! COPYRIGHT
!! Copyright (C) 2001-2018 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_multibinit_unittest
  use defs_basis
  use m_abicore
  use m_errors

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_abstract_potential, only: abstract_potential_t

  implicit none
  !!***

contains

  subroutine mb_test_main()
    call mb_test1()
  end subroutine mb_test_main
  !-------------------------------------------------

  subroutine mb_test1()
    use m_supercell_maker, only: scmaker_unittest
    call scmaker_unittest()
  end subroutine mb_test1

  !-------------------------------------------------
end module m_multibinit_unittest
